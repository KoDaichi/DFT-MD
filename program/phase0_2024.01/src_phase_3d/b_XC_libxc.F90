#ifdef __EDA__
subroutine ex_lda_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, &
     &                   wos, exc, dFx_drho, exc_on_a_grid_wk, ist, ien )
#else
subroutine ex_lda_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, &
     &                   wos, exc, dFx_drho, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI

  use m_Control_Parameters,  only : xc_func_exch
  use xc_f03_lib_m

#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer(8) :: npoint
  integer :: i, j, is
  integer :: istart, iend
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: exc0(:)

  istart = ista_r;   iend = iend_r
  istart = ist;      iend   = ien

  npoint = iend -istart +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( vrho(nspin,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +istart -1
        rho(1,i) = chgrhr_l(j,1)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
        End Do
     else
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)
        End Do
     endif
  End if

  call xc_f03_lda_exc_vxc( xc_func_exch, npoint, rho, exc0, vrho )
!
  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=istart, iend
        j = i -istart +1
        dFx_drho(i,1) = vrho(1,j)
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=istart, iend
           j = i -istart +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(1,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=istart, iend
           j = i -istart +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(2,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

#ifdef __EDA__
  if(sw_eda==ON) then
     Do is=1, ispin
        Do i=ista_r, iend_r
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) +exc0(i) *chgrhr_l(i,is)
        End Do
     End Do
  endif
#endif

  deallocate( rho );   deallocate( vrho );  deallocate( exc0 )

end subroutine ex_lda_libxc

#ifdef __EDA__
subroutine cr_lda_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, &
     &                   wos, exc, dF_drho, exc_on_a_grid_wk, ist, ien )
#else
subroutine cr_lda_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, &
     &                   wos, exc, dF_drho, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_corr
  use xc_f03_lib_m
  
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(out) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer(8) :: npoint
  integer :: i, j, is
  integer :: istart, iend
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: exc0(:)

  istart = ista_r;   iend = iend_r
  istart = ist;      iend   = ien

  npoint = iend -istart +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( vrho(nspin,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +istart -1
        rho(1,i) = chgrhr_l(j,1)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
        End Do
     else
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1)
           rho(2,i) = chgrhr_l(j,2)
        End Do
     endif
  End if

  call xc_f03_lda_exc_vxc( xc_func_corr, npoint, rho, exc0, vrho )

!!  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=istart, iend
        j = i -istart +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=istart, iend
           j = i -istart +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(1,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=istart, iend
           j = i -istart +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(2,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

#ifdef __EDA__
  if(sw_eda==ON) then
     Do is=1, ispin
        Do i=ista_r, iend_r
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0(i) *chgrhr_l(i,is)
        End Do
     End Do
  endif
#endif

  deallocate( rho );  deallocate( vrho );   deallocate( exc0 )

end subroutine cr_lda_libxc

#ifdef __EDA__
subroutine ex_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
     &                   wos, exc, dFx_drho, dFx_dgradrho, exc_on_a_grid_wk, ist, ien )
#else
subroutine ex_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
     &                   wos, exc, dFx_drho, dFx_dgradrho, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI

  use m_Control_Parameters,  only : xc_func_exch
  use xc_f03_lib_m
  use m_Parallelization,    only : mype

#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
!  integer,intent(in), optional :: ist, ien
  integer,intent(in) :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer(8) :: npoint
  integer :: i, j, is
  integer :: istart, iend
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: exc0(:)

  istart = ista_r;   iend = iend_r
  istart = ist;      iend   = ien

  npoint = iend -istart +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( sigma(2*nspin-1,npoint) )
  allocate( vrho(nspin,npoint) )
  allocate( vsigma(2*nspin-1,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +istart -1
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
           sigma(1,i) = grad_rho(j,1)**2 /4.0d0
           sigma(2,i) = sigma(1,j)
           sigma(3,i) = sigma(1,j)
        End Do
     else
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)

           c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
           sigma(1,i) = c1
           sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
           sigma(3,i) = c2
        End Do
     endif
  End if

  call xc_f03_gga_exc_vxc( xc_func_exch, npoint, rho, sigma, exc0, vrho, vsigma )
!
  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=istart, iend
        j = i -istart +1
        dFx_drho(i,1) = vrho(1,j)
        dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=istart, iend
           j = i -istart +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(1,j)
           dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
           dFx_dgradrho(i,2) = vsigma(1,j) *2.0d0
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=istart, iend
           j = i -istart +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(2,j)
           dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
           dFx_dgradrho(i,2) = vsigma(3,j) *2.0d0
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

#ifdef __EDA__
  if(sw_eda==ON) then
     Do is=1, ispin
        Do i=ista_r, iend_r
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) +exc0(i) *chgrhr_l(i,is)
        End Do
     End Do
  endif
#endif

  deallocate( rho );  deallocate( sigma )
  deallocate( vrho ); deallocate( vsigma )
  deallocate( exc0 )

end subroutine ex_gga_libxc

#ifdef __EDA__
subroutine cr_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
     &                   wos, exc, dF_drho, dF_dgradrho, exc_on_a_grid_wk, ist, ien )
#else
subroutine cr_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
     &                   wos, exc, dF_drho, dF_dgradrho, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_corr
  use xc_f03_lib_m
  
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(out) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_r:iend_r,1)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer(8) :: npoint
  integer :: i, j, is
  integer :: istart, iend
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: exc0(:)

  istart = ista_r;   iend = iend_r
  istart = ist;      iend   = ien

  npoint = iend -istart +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( sigma(2*nspin-1,npoint) )
  allocate( vrho(nspin,npoint) )
  allocate( vsigma(2*nspin-1,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +istart -1
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
           sigma(1,i) = grad_rho(j,1)**2 /4.0d0
           sigma(2,i) = sigma(1,i)
           sigma(3,i) = sigma(1,i)
        End Do
     else
        Do i=1, npoint
           j = i +istart -1
           rho(1,i) = chgrhr_l(j,1)
           rho(2,i) = chgrhr_l(j,2)

           c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
           sigma(1,i) = c1
           sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
           sigma(3,i) = c2
        End Do
     endif
  End if

  call xc_f03_gga_exc_vxc( xc_func_corr, npoint, rho, sigma, exc0, vrho, vsigma )

!!  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=istart, iend
        j = i -istart +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        dF_dgradrho(i,1) = vsigma(1,j) *2.0d0
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=istart, iend
           j = i -istart +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(1,j)
           dF_dgradrho(i,1) = vsigma(1,j) *2.0d0
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=istart, iend
           j = i -istart +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(2,j)
           dF_dgradrho(i,1) = ( vsigma(1,j) +vsigma(3,j) )    ! ???
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

#ifdef __EDA__
  if(sw_eda==ON) then
     Do is=1, ispin
        Do i=ista_r, iend_r
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0(i) *chgrhr_l(i,is)
        End Do
     End Do
  endif
#endif

  deallocate( rho );  deallocate( sigma )
  deallocate( vrho ); deallocate( vsigma )
  deallocate( exc0 )

end subroutine cr_gga_libxc

subroutine ex_mgga_libxc( nspin, ispin, ista_r, iend_r, &
     &                    chgrhr_l, grad_rho, grad_trho, lapl_rho, ekin_dens, &
     &                    wos, exc,  &
     &                    dFx_drho, dFx_dgradrho, dFx_dlaplrho, dFx_dtau )
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_exch, xc_info_exch, xc_flag_exch, &
       &                            xc_name_exch, val_c_tb09
  use xc_f03_lib_m

  implicit none

  integer, intent(in) :: nspin, ispin, ista_r, iend_r
  real(kind=DP), intent(in) :: chgrhr_l( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: grad_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: grad_trho( ista_r:iend_r )
  real(kind=DP), intent(in) :: lapl_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: ekin_dens( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: wos(ista_r:iend_r)
  real(kind=DP), intent(out) :: dFx_drho( ista_r:iend_r, nspin )
  real(kind=DP), intent(out) :: dFx_dgradrho( ista_r:iend_r, nspin )
  real(kind=DP), intent(out) :: dFx_dlaplrho( ista_r:iend_r, nspin )
  real(kind=DP), intent(out) :: dFx_dtau( ista_r:iend_r, nspin )
  real(kind=DP),intent(inout) :: exc
!
  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: lapla(:,:), tau(:,:), vlapla(:,:), vtau(:,:)
  real(kind=DP), allocatable :: exc0(:)

  npoint = iend_r -ista_r +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( sigma(2*nspin-1,npoint) )
  allocate( lapla(nspin,npoint) )
  allocate( tau(nspin,npoint) )

  allocate( vrho(nspin,npoint) )
  allocate( vsigma(2*nspin-1,npoint) )
  allocate( vlapla(nspin,npoint) )
  allocate( vtau(nspin,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +ista_r -1
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
        lapla(1,i) = lapl_rho(j,1)
        tau(1,i)   = ekin_dens(j,1)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +ista_r -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
           sigma(1,i) = grad_rho(j,1)**2 /4.0d0
           sigma(2,i) = sigma(1,i)
           sigma(3,i) = sigma(1,i)
           lapla(1,i) = lapl_rho(j,1) /2.0d0
           lapla(2,i) = lapla(1,i)
           tau(1,i)   = ekin_dens(j,1) /2.0d0
           tau(2,i)   = tau(1,i)
        End Do
     else
        Do i=1, npoint
           j = i +ista_r -1
           rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)

           c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
           sigma(1,i) = c1
           sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
           sigma(3,i) = c2

           lapla(1,i) = lapl_rho(j,1);    lapla(2,i) = lapl_rho(j,2) 
           tau(1,i)   = ekin_dens(j,1);   tau(2,i)   = ekin_dens(j,2)
        End Do
     endif
  End if

  if ( xc_name_exch == "mgga_x_tb09" ) then
     call xc_f03_func_set_ext_params( xc_func_exch, [val_c_tb09] )
  endif

  if ( xc_flag_exch(0) == 0 ) then         ! pot_only
     call xc_f03_mgga_vxc( xc_func_exch, npoint, rho, sigma, lapla, tau, &
          &                vrho, vsigma, vlapla, vtau )
     exc0 = 0.0d0
  else
     call xc_f03_mgga_exc_vxc( xc_func_exch, npoint, rho, sigma, lapla, tau, &
          &                    exc0, vrho, vsigma, vlapla, vtau )
  endif
!
  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ista_r, iend_r
        j = i -ista_r +1
        dFx_drho(i,1) = vrho(1,j)
        dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
        dFx_dlaplrho(i,1) = vlapla(1,j)
        dFx_dtau(i,1) = vtau(1,j)
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=ista_r, iend_r
           j = i -ista_r +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(1,j)
           dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
           dFx_dgradrho(i,2) = vsigma(1,j) *2.0d0
           dFx_dlaplrho(i,1) = vlapla(1,j)
           dFx_dlaplrho(i,2) = vlapla(1,j)
           dFx_dtau(i,1) = vtau(1,j)
           dFx_dtau(i,2) = vtau(1,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=ista_r, iend_r
           j = i -ista_r +1
           dFx_drho(i,1) = vrho(1,j)
           dFx_drho(i,2) = vrho(2,j)
           dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0
           dFx_dgradrho(i,2) = vsigma(3,j) *2.0d0
           dFx_dlaplrho(i,1) = vlapla(1,j)
           dFx_dlaplrho(i,2) = vlapla(2,j)
           dFx_dtau(i,1) = vtau(1,j)
           dFx_dtau(i,2) = vtau(2,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

  deallocate( rho );  deallocate( sigma );  deallocate( lapla );  deallocate( tau )
  deallocate( vrho ); deallocate( vsigma ); deallocate( vlapla ); deallocate( vtau )
  deallocate( exc0 )
!
end subroutine ex_mgga_libxc

subroutine cr_mgga_libxc( nspin, ispin, ista_r, iend_r, &
     &                    chgrhr_l, grad_rho, grad_trho, lapl_rho, ekin_dens, &
     &                    wos, exc,  &
     &                    dF_drho, dF_dgradrho, dF_dlaplrho, dF_dtau )
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_corr
  use xc_f03_lib_m

  implicit none

  integer, intent(in) :: nspin, ispin, ista_r, iend_r
  real(kind=DP), intent(in) :: chgrhr_l( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: grad_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: grad_trho( ista_r:iend_r )
  real(kind=DP), intent(in) :: lapl_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: ekin_dens( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: wos(ista_r:iend_r)
  real(kind=DP), intent(inout) :: dF_drho( ista_r:iend_r, nspin )
  real(kind=DP), intent(inout) :: dF_dgradrho( ista_r:iend_r, 1 )
  real(kind=DP), intent(inout) :: dF_dlaplrho( ista_r:iend_r, nspin )
  real(kind=DP), intent(inout) :: dF_dtau( ista_r:iend_r, nspin )
  real(kind=DP),intent(inout) :: exc
!
  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: lapla(:,:), tau(:,:), vlapla(:,:), vtau(:,:)
  real(kind=DP), allocatable :: exc0(:)

!
  npoint = iend_r -ista_r +1
  if ( npoint < 1 ) return

  allocate( rho(nspin,npoint) )
  allocate( sigma(2*nspin-1,npoint) )
  allocate( lapla(nspin,npoint) )
  allocate( tau(nspin,npoint) )

  allocate( vrho(nspin,npoint) )
  allocate( vsigma(2*nspin-1,npoint) )
  allocate( vlapla(nspin,npoint) )
  allocate( vtau(nspin,npoint) )
  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = i +ista_r -1
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
        lapla(1,i) = lapl_rho(j,1)
        tau(1,i)   = ekin_dens(j,1)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=1, npoint
           j = i +ista_r -1
           rho(1,i) = chgrhr_l(j,1) /2.0d0
           rho(2,i) = rho(1,i)
           sigma(1,i) = grad_rho(j,1)**2 /4.0d0
           sigma(2,i) = sigma(1,i)
           sigma(3,i) = sigma(1,i)
           lapla(1,i) = lapl_rho(j,1) /2.0d0
           lapla(2,i) = lapla(1,i)
           tau(1,i)   = ekin_dens(j,1) /2.0d0
           tau(2,i)   = tau(1,i)
        End Do
     else
        Do i=1, npoint
           j = i +ista_r -1
           rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)

           c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
           sigma(1,i) = c1
           sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
           sigma(3,i) = c2

           lapla(1,i) = lapl_rho(j,1);    lapla(2,i) = lapl_rho(j,2) 
           tau(1,i)   = ekin_dens(j,1);   tau(2,i)   = ekin_dens(j,2)
        End Do
     endif
  End if

  call xc_f03_mgga_exc_vxc( xc_func_corr, npoint, rho, sigma, lapla, tau, &
       &                    exc0, vrho, vsigma, vlapla, vtau )

!!  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ista_r, iend_r
        j = i -ista_r +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        dF_dgradrho(i,1) = vsigma(1,j) *2.0d0
        dF_dlaplrho(i,1) = dF_dlaplrho(i,1) +vlapla(1,j)
        dF_dtau(i,1)     = dF_dtau(i,1)     +vtau(1,j)
        exc = exc +wos(i) *exc0(j) *rho(1,j)
     End Do
  else
     if ( ispin == 1 ) then
        Do i=ista_r, iend_r
           j = i -ista_r +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(1,j)
           dF_dgradrho(i,1) = vsigma(1,j) *2.0d0
           dF_dlaplrho(i,1) = dF_dlaplrho(i,1) +vlapla(1,j)
           dF_dlaplrho(i,2) = dF_dlaplrho(i,2) +vlapla(1,j)
           dF_dtau(i,1)     = dF_dtau(i,1)     +vtau(1,j)
           dF_dtau(i,2)     = dF_dtau(i,2)     +vtau(1,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j)*2.0d0 )
        End Do
     else
        Do i=ista_r, iend_r
           j = i -ista_r +1
           dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
           dF_drho(i,2) = dF_drho(i,2) +vrho(2,j)
           dF_dgradrho(i,1) = vsigma(1,j) +vsigma(3,j)

           dF_dlaplrho(i,1) = dF_dlaplrho(i,1) +vlapla(1,j)
           dF_dlaplrho(i,2) = dF_dlaplrho(i,2) +vlapla(2,j)
           dF_dtau(i,1)     = dF_dtau(i,1) +vtau(1,j)
           dF_dtau(i,2)     = dF_dtau(i,2) +vtau(2,j)
           exc = exc +wos(i) *exc0(j) *( rho(1,j) +rho(2,j) )
        End Do
     endif
  End if

  deallocate( rho );  deallocate( sigma );  deallocate( lapla );  deallocate( tau )
  deallocate( vrho ); deallocate( vsigma ); deallocate( vlapla ); deallocate( vtau )
  deallocate( exc0 )

end subroutine cr_mgga_libxc

! ======== PAW Surface Harmonics ====== 
subroutine ex_gga_paw_libxc( nrc, dnr, nspin, chgrhr_l, grad_rho, grad_trho, exc, &
     &                       dFx_drho, dFx_dgradrho, &
     &                       dFx_drr,  dFx_drg,  dFx_dgg, &
     &                       dFx_drrr, dFx_drrg, dFx_drgg, dFx_dggg, &
     &                       ista_nrc, iend_nrc, ist, ien )
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_exch, xc_flag_exch
  use m_Parallelization,   only : mype
  use xc_f03_lib_m

  implicit none

  integer,intent(in)        :: nrc, dnr, nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in)  :: grad_trho(ista_nrc:iend_nrc)

  real(kind=DP),intent(out) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFx_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(ista_nrc:iend_nrc,nspin)
  integer :: ista_nrc, iend_nrc, ist, ien

  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: v2rho2(:,:), v2rhosigma(:,:), v2sigma2(:,:)
  real(kind=DP), allocatable :: v3rho3(:,:), v3rho2sigma(:,:), v3rhosigma2(:,:), &
       &                        v3sigma3(:,:)
  real(kind=DP), allocatable :: v4rho4(:,:), v4rho3sigma(:,:), v4rho2sigma2(:,:), &
       &                        v4rhosigma3(:,:), v4sigma4(:,:)
  real(kind=DP), allocatable :: exc0(:)

  npoint = (ien -ist +1)/dnr
  if ( npoint < 1 ) return

  if ( nspin == 1 ) then
     allocate( rho(1,npoint) );      allocate( sigma(1,npoint) )
     allocate( vrho(1,npoint) );     allocate( vsigma(1,npoint) )
     allocate( v2rho2(1,npoint) );   allocate( v2rhosigma(1,npoint) )
     allocate( v2sigma2(1,npoint) )
#ifdef KXC_LXC_AVAILABLE
     if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
        allocate( v3rho3(1,npoint) );   allocate( v3rho2sigma(1,npoint) )
        allocate( v3rhosigma2(1,npoint) );  
        allocate( v3sigma3(1,npoint) )
        allocate( v4rho4(1,npoint) );   allocate( v4rho3sigma(1,npoint) )
        allocate( v4rho2sigma2(1,npoint) );  
        allocate( v4rhosigma3(1,npoint) );
        allocate( v4sigma4(1,npoint) )
     endif
#endif
  else
     allocate( rho(2,npoint) );      allocate( sigma(3,npoint) )
     allocate( vrho(2,npoint) );     allocate( vsigma(3,npoint) )
     allocate( v2rho2(3,npoint) );   allocate( v2rhosigma(6,npoint) )
     allocate( v2sigma2(6,npoint) )
#ifdef KXC_LXC_AVAILABLE
     if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
        allocate( v3rho3(4,npoint) );   allocate( v3rho2sigma(9,npoint) )
        allocate( v3rhosigma2(12,npoint) );  
        allocate( v3sigma3(10,npoint) )
        allocate( v4rho4(5,npoint) );   allocate( v4rho3sigma(12,npoint) )
        allocate( v4rho2sigma2(15,npoint) );  
        allocate( v4rhosigma3(20,npoint) );
        allocate( v4sigma4(15,npoint) )
     endif
#endif
  endif

  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
     End Do
  else
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)
        
        c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
        sigma(1,i) = c1
        sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
        sigma(3,i) = c2
     End Do
  End if

#ifdef KXC_LXC_AVAILABLE
  if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
     call xc_f03_gga( xc_func_exch, npoint, rho, sigma, exc0, vrho, vsigma, &
          &           v2rho2, v2rhosigma,  v2sigma2, &
          &           v3rho3, v3rho2sigma, v3rhosigma2,  v3sigma3, &
          &           v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 )
  else
     call xc_f03_gga_exc_vxc_fxc( xc_func_exch, npoint, rho, sigma, exc0, &
          &           vrho, vsigma, &
          &           v2rho2, v2rhosigma, v2sigma2 )
  endif
#else
  call xc_f03_gga_exc_vxc_fxc( xc_func_exch, npoint, rho, sigma, exc0, &
       &           vrho, vsigma, &
       &           v2rho2, v2rhosigma, v2sigma2 )
#endif

  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dFx_drho(i,1) = vrho(1,j)
        dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0 *grad_rho(i,1)
        exc(i) = exc(i) +exc0(j) *rho(1,j)

        dFx_drr(i,1) = v2rho2(1,j)
        dFx_drg(i,1) = v2rhosigma(1,j) *2.0d0 *grad_rho(i,1)
        dFx_dgg(i,1) = vsigma(1,j) *2.0d0 &
             &         +v2sigma2(1,j) *4.0d0 *grad_rho(i,1)**2 

#ifdef KXC_LXC_AVAILABLE
        if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
           dFx_drrr(i,1) = v3rho3(1,j)
           dFx_drrg(i,1) = v3rho2sigma(1,j) *2.0d0 *grad_rho(i,1)
           dFx_drgg(i,1) = v2rhosigma(1,j) *2.0d0 &
                &          +v3rhosigma2(1,j) *4.0d0 *grad_rho(i,1)**2
           dFx_dggg(i,1) = vsigma2(1,j) *12.0d0 *grad_rho(i,1) &
                &          +v3sigma3(1,j) *8.0d0 *grad_rho(i,1)**3
        endif
#endif
     End Do
  else
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dFx_drho(i,1) = vrho(1,j)
        dFx_drho(i,2) = vrho(2,j)
        dFx_dgradrho(i,1) = vsigma(1,j) *2.0d0 *grad_rho(i,1)
        dFx_dgradrho(i,2) = vsigma(3,j) *2.0d0 *grad_rho(i,2)

        exc(i) = exc(i) +exc0(j) *( rho(1,j) +rho(2,j) )
        dFx_drr(i,1) = v2rho2(1,j)
        dFx_drr(i,2) = v2rho2(3,j)
        dFx_drg(i,1) = v2rhosigma(1,j) *2.0d0 *grad_rho(i,1)
        dFx_drg(i,2) = v2rhosigma(6,j) *2.0d0 *grad_rho(i,2)
        dFx_dgg(i,1) = vsigma(1,j) *2.0d0 &
             &         +v2sigma2(1,j) *4.0d0 *grad_rho(i,1)**2 
        dFx_dgg(i,2) = vsigma(3,j) *2.0d0 &
             &         +v2sigma2(6,j) *4.0d0 *grad_rho(i,2)**2 

#ifdef KXC_LXC_AVAILABLE
        if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
           dFx_drrr(i,1) = v3rho3(1,j)
           dFx_drrr(i,2) = v3rho3(4,j)
           dFx_drrg(i,1) = v3rho2sigma(1,j) *2.0d0 *grad_rho(i,1)
           dFx_drrg(i,2) = v3rho2sigma(9,j) *2.0d0 *grad_rho(i,2)
           
           dFx_drgg(i,1) = v2rhosigma(1,j) *2.0d0 &
                &          +v3rhosigma2(1,j) *4.0d0 *grad_rho(i,1)**2
           dFx_drgg(i,2) = v2rhosigma(6,j) *2.0d0 &
                &          +v3rhosigma2(12,j) *4.0d0 *grad_rho(i,2)**2
           dFx_dggg(i,1) = v2sigma2(1,j) *12.0d0 *grad_rho(i,1) &
                &          +v3sigma3(1,j) *8.0d0 *grad_rho(i,1)**3
           dFx_dggg(i,2) = v2sigma2(6,j) *12.0d0 *grad_rho(i,2) &
                &          +v3sigma3(10,j) *8.0d0 *grad_rho(i,2)**3
        endif
#endif
     End Do
  End if

  deallocate( rho );  deallocate( sigma )
  deallocate( vrho ); deallocate( vsigma )
  deallocate( v2rho2 ); deallocate( v2rhosigma );  deallocate( v2sigma2 )

#ifdef KXC_LXC_AVAILABLE
  if ( xc_flag_exch(3) == 1 .and. xc_flag_exch(4) == 1 ) then
     deallocate( v3rho3 ); deallocate( v3rho2sigma );  deallocate( v3rhosigma2 );
     deallocate( v3sigma3 )
     deallocate( v4rho4 ); deallocate( v4rho3sigma );  deallocate( v4rho2sigma2 );
     deallocate( v4rhosigma3 ); deallocate( v4sigma4 )
  endif
#endif

  deallocate( exc0 )

end subroutine ex_gga_paw_libxc

subroutine cr_gga_paw_libxc( nrc, dnr, nspin, chgrhr_l, grad_rho, grad_trho, exc, &
     &                       dF_drho,  dF_dgradrho, dFc_daa,  dFc_dbb,  dFc_dgg ,&
     &                       dFc_dab,  dFc_dag,     dFc_dbg,  dFc_daaa, dFc_dbbb, &
     &                       dFc_dggg, dFc_daab,    dFc_daag, dFc_dabb, dFc_dbbg, &
     &                       dFc_dagg, dFc_dbgg, dFc_dabg, &
     &                       ista_nrc, iend_nrc, ist, ien )
  use m_Const_Parameters,  only : PAI,DP
  use m_Control_Parameters,  only : xc_func_corr, xc_flag_corr
  use xc_f03_lib_m

  implicit none

  integer,intent(in)        :: nrc, dnr, nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in) :: grad_rho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(ista_nrc:iend_nrc)
  real(kind=DP),intent(inout) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dF_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daaa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dggg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dagg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabg(ista_nrc:iend_nrc)
  integer :: ista_nrc, iend_nrc, ist, ien
  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), sigma(:,:), vrho(:,:), vsigma(:,:)
  real(kind=DP), allocatable :: v2rho2(:,:), v2rhosigma(:,:), v2sigma2(:,:)
  real(kind=DP), allocatable :: v3rho3(:,:), v3rho2sigma(:,:), v3rhosigma2(:,:), &
       &                        v3sigma3(:,:)
  real(kind=DP), allocatable :: v4rho4(:,:), v4rho3sigma(:,:), v4rho2sigma2(:,:), &
       &                        v4rhosigma3(:,:), v4sigma4(:,:)
  real(kind=DP), allocatable :: exc0(:)

  npoint = (ien -ist +1)/dnr
  if ( npoint < 1 ) return

  if ( nspin == 1 ) then
     allocate( rho(1,npoint) );      allocate( sigma(1,npoint) )
     allocate( vrho(1,npoint) );     allocate( vsigma(1,npoint) )
     allocate( v2rho2(1,npoint) );   allocate( v2rhosigma(1,npoint) )
     allocate( v2sigma2(1,npoint) )
#ifdef KXC_LXC_AVAILABLE
     if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
        allocate( v3rho3(1,npoint) );   allocate( v3rho2sigma(1,npoint) )
        allocate( v3rhosigma2(1,npoint) );  
        allocate( v3sigma3(1,npoint) )
        allocate( v4rho4(1,npoint) );   allocate( v4rho3sigma(1,npoint) )
        allocate( v4rho2sigma2(1,npoint) );  
        allocate( v4rhosigma3(1,npoint) );
        allocate( v4sigma4(1,npoint) )
     endif
#endif
  else
     allocate( rho(2,npoint) );      allocate( sigma(3,npoint) )
     allocate( vrho(2,npoint) );     allocate( vsigma(3,npoint) )
     allocate( v2rho2(3,npoint) );   allocate( v2rhosigma(6,npoint) )
     allocate( v2sigma2(6,npoint) )
#ifdef KXC_LXC_AVAILABLE
     if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
        allocate( v3rho3(4,npoint) );   allocate( v3rho2sigma(9,npoint) )
        allocate( v3rhosigma2(12,npoint) );  
        allocate( v3sigma3(10,npoint) )
        allocate( v4rho4(5,npoint) );   allocate( v4rho3sigma(12,npoint) )
        allocate( v4rho2sigma2(15,npoint) );  
        allocate( v4rhosigma3(20,npoint) );
        allocate( v4sigma4(15,npoint) )
     endif
#endif
  endif

  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1)
        sigma(1,i) = grad_rho(j,1)**2
     End Do
  else
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)
        
        c1 = grad_rho(j,1)**2;         c2 = grad_rho(j,2)**2
        sigma(1,i) = c1
        sigma(2,i) = ( grad_trho(j)**2 -c1 -c2 ) /2.0d0
        sigma(3,i) = c2
     End Do
  End if

#ifdef KXC_LXC_AVAILABLE
  if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
     call xc_f03_gga( xc_func_corr, npoint, rho, sigma, exc0, vrho, vsigma, &
          &           v2rho2, v2rhosigma,  v2sigma2, &
          &           v3rho3, v3rho2sigma, v3rhosigma2,  v3sigma3, &
          &           v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 )
  else
     call xc_f03_gga_exc_vxc_fxc( xc_func_corr, npoint, rho, sigma, exc0, &
          &           vrho, vsigma, &
          &           v2rho2, v2rhosigma,  v2sigma2 )
  endif
#else
  call xc_f03_gga_exc_vxc_fxc( xc_func_corr, npoint, rho, sigma, exc0, &
       &           vrho, vsigma, &
       &           v2rho2, v2rhosigma,  v2sigma2 )
#endif
!
!  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        dF_dgradrho(i) = vsigma(1,j) *2.0d0 *grad_rho(i,1)
        exc(i) = exc(i) +exc0(j) *rho(1,j)

        dFc_daa(i) = v2rho2(1,j)
        dFc_dag(i) = v2rhosigma(1,j) *2.0d0 *grad_rho(i,1)
        dFc_dgg(i) = vsigma(1,j) *2.0d0 &
             &       +v2sigma2(1,j) *4.0d0 *grad_rho(i,1)**2

#ifdef KXC_LXC_AVAILABLE
        if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
           dFc_daaa(i) = v3rho3(1,j)
           dFc_daag(i) = v3rho2sigma(1,j) *2.0d0 *grad_rho(i,1)
           dFc_dagg(i) = v2rhosigma(1,j) *2.0d0 &
                &        +v3rhosigma2(1,j) *4.0d0 *grad_rho(i,1)**2
           dFc_dggg(i) = vsigma2(1,j) *12.0d0 *grad_rho(i,1) &
                &        +v3sigma3(1,j) *8.0d0 *grad_rho(i,1)**3
        endif
#endif
     End Do
  else
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        dF_drho(i,2) = dF_drho(i,2) +vrho(2,j)

        dF_dgradrho(i) = ( vsigma(1,j) +vsigma(3,j) ) *grad_trho(i)

        exc(i) = exc(i) +exc0(j) *( rho(1,j) +rho(2,j) )

        dFc_daa(i) = v2rho2(1,j)
        dFc_dab(i) = v2rho2(2,j)
        dFc_dbb(i) = v2rho2(3,j)
        dFc_dag(i) = v2rhosigma(1,j) *2.0d0 *( grad_trho(i) )
        dFc_dbg(i) = v2rhosigma(6,j) *2.0d0 *( grad_trho(i) )

! is this correct ??
        dFc_dgg(i) = ( vsigma(1,j) +vsigma(3,j) )&  
             &       +v2sigma2(1,j) *4.0d0 *grad_rho(i,1) *grad_trho(i) &
             &       +v2sigma2(6,j) *4.0d0 *grad_rho(i,2) *grad_trho(i)

#ifdef KXC_LXC_AVAILABLE
        if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
           stop "KXC_LXC are not supported"
           dFc_daaa(i) = v3rho3(1,j)
           dFc_daab(i) = v3rho3(2,j)
           dFc_dabb(i) = v3rho3(3,j)
           dFc_dbbb(i) = v3rho3(4,j)
           
           dFc_daag(i) = v3rho2sigma(1,j) *2.0d0 *grad_rho(i,1)
           dFc_dagg(i) = v2rhosigma(1,j) *2.0d0 &
                &        +v3rhosigma2(1,j) *4.0d0 *grad_rho(i,1)**2
           dFc_dggg(i) = vsigma2(1,j) *12.0d0 *grad_rho(i,1) &
                &        +v3sigma3(1,j) *8.0d0 *grad_rho(i,1)**3
        endif
#endif
     End Do
  End if

  deallocate( rho );  deallocate( sigma )
  deallocate( vrho ); deallocate( vsigma )
  deallocate( v2rho2 ); deallocate( v2rhosigma );  deallocate( v2sigma2 )
#ifdef KXC_LXC_AVAILABLE
  if ( xc_flag_corr(3) == 1 .and. xc_flag_corr(4) == 1 ) then
     deallocate( v3rho3 ); deallocate( v3rho2sigma );  deallocate( v3rhosigma2 );
     deallocate( v3sigma3 )
     deallocate( v4rho4 ); deallocate( v4rho3sigma );  deallocate( v4rho2sigma2 );
     deallocate( v4rhosigma3 ); deallocate( v4sigma4 )
  endif
#endif
  deallocate( exc0 )

end subroutine cr_gga_paw_libxc

subroutine ex_lda_paw_libxc( nrc, dnr, nspin, chgrhr_l, exc, &
     &                       dFx_drho, dFx_drr, dFx_drrr, &
     &                       ista_nrc, iend_nrc, ist, ien )
  use m_Const_Parameters,  only : DP,PAI
  use m_Control_Parameters,  only : xc_func_exch
  use xc_f03_lib_m

  implicit none

  integer,intent(in)        :: nrc, dnr, nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)

  real(kind=DP),intent(out) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFx_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(ista_nrc:iend_nrc,nspin)
  integer :: ista_nrc, iend_nrc, ist, ien

  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:), vrho(:,:)
  real(kind=DP), allocatable :: v2rho2(:,:)
  real(kind=DP), allocatable :: v3rho3(:,:)
  real(kind=DP), allocatable :: v4rho4(:,:)
  real(kind=DP), allocatable :: exc0(:)

  npoint = (ien -ist +1)/dnr
  if ( npoint < 1 ) return

  if ( nspin == 1 ) then
     allocate( rho(1,npoint) ); 
     allocate( vrho(1,npoint) );
     allocate( v2rho2(1,npoint) );
#ifdef KXC_LXC_AVAILABLE
     allocate( v3rho3(1,npoint) );
     allocate( v4rho4(1,npoint) );
#endif
  else
     allocate( rho(2,npoint) );
     allocate( vrho(2,npoint) );
     allocate( v2rho2(3,npoint) );
#ifdef KXC_LXC_AVAILABLE
     allocate( v3rho3(4,npoint) );
     allocate( v4rho4(5,npoint) );
#endif
  endif

  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1)
     End Do
  else
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)
     End Do
  End if

#ifdef KXC_LXC_AVAILABLE
  call xc_f03_lda( xc_func_exch, npoint, rho, exc0, vrho, &
       &           v2rho2, v3rho3, v4rho4 )
#else
  call xc_f03_lda_exc_vxc_fxc( xc_func_exch, npoint, rho, exc0, &
       &           vrho, v2rho2 )
#endif

  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dFx_drho(i,1) = vrho(1,j)
        exc(i) = exc(i) +exc0(j) *rho(1,j)

        dFx_drr(i,1) = v2rho2(1,j)
#ifdef KXC_LXC_AVAILABLE
        dFx_drrr(i,1) = v3rho3(1,j)
#endif
     End Do
  else
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dFx_drho(i,1) = vrho(1,j)
        dFx_drho(i,2) = vrho(2,j)

        exc(i) = exc(i) +exc0(j) *( rho(1,j) +rho(2,j) )
        dFx_drr(i,1) = v2rho2(1,j)
        dFx_drr(i,2) = v2rho2(3,j)
#ifdef KXC_LXC_AVAILABLE
        dFx_drrr(i,1) = v3rho3(1,j)
        dFx_drrr(i,2) = v3rho3(4,j)
#endif
     End Do
  End if

  deallocate( rho );
  deallocate( vrho )
  deallocate( v2rho2 )

#ifdef KXC_LXC_AVAILABLE
  deallocate( v3rho3 )
  deallocate( v4rho4 )
#endif

  deallocate( exc0 )
  
end subroutine ex_lda_paw_libxc

subroutine cr_lda_paw_libxc( nrc, dnr, nspin, chgrhr_l, exc, &
     &                       dF_drho,  dFc_daa,  dFc_dbb, &
     &                       dFc_dab,  dFc_daaa, dFc_dbbb, &
     &                       dFc_daab, dFc_dabb, &
     &                       ista_nrc, iend_nrc, ist, ien )
  use m_Const_Parameters,  only : PAI,DP
  use m_Control_Parameters,  only : xc_func_corr
  use xc_f03_lib_m

  implicit none

  integer,intent(in)        :: nrc, dnr, nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(inout) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dF_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFc_daa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daaa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabb(ista_nrc:iend_nrc)
  integer :: ista_nrc, iend_nrc, ist, ien
  integer(8) :: npoint
  integer :: i, j, is
  real(kind=DP) :: c1, c2
  real(kind=DP), allocatable :: rho(:,:)
  real(kind=DP), allocatable :: vrho(:,:)
  real(kind=DP), allocatable :: v2rho2(:,:)
  real(kind=DP), allocatable :: v3rho3(:,:)
  real(kind=DP), allocatable :: v4rho4(:,:)
  real(kind=DP), allocatable :: exc0(:)

  npoint = (ien -ist +1)/dnr
  if ( npoint < 1 ) return

  if ( nspin == 1 ) then
     allocate( rho(1,npoint) )
     allocate( vrho(1,npoint) )
     allocate( v2rho2(1,npoint) )
#ifdef KXC_LXC_AVAILABLE
     allocate( v3rho3(1,npoint) )
     allocate( v4rho4(1,npoint) ); 
#endif
  else
     allocate( rho(2,npoint) );    
     allocate( vrho(2,npoint) );   
     allocate( v2rho2(3,npoint) ); 
#ifdef KXC_LXC_AVAILABLE
     allocate( v3rho3(4,npoint) ); 
     allocate( v4rho4(5,npoint) );
#endif
  endif

  allocate( exc0(npoint) )

  if ( nspin == 1 ) then
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1)
     End Do
  else
     Do i=1, npoint
        j = ist +(i-1) *dnr
        rho(1,i) = chgrhr_l(j,1);      rho(2,i) = chgrhr_l(j,2)
     End Do
  End if

#ifdef KXC_LXC_AVAILABLE
  call xc_f03_lda( xc_func_corr, npoint, rho, exc0, vrho, &
       &           v2rho2, v3rho3, v4rho4 )
#else
  call xc_f03_lda_exc_vxc_fxc( xc_func_corr, npoint, rho, exc0, &
       &           vrho, v2rho2 )
#endif
!
!  exc = 0.0d0
  if ( nspin == 1 ) then
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        exc(i) = exc(i) +exc0(j) *rho(1,j)

        dFc_daa(i) = v2rho2(1,j)
#ifdef KXC_LXC_AVAILABLE
        dFc_daaa(i) = v3rho3(1,j)
#endif
     End Do
  else
     Do i=ist, ien, dnr
        j = ( i-ist ) /dnr +1
        dF_drho(i,1) = dF_drho(i,1) +vrho(1,j)
        dF_drho(i,2) = dF_drho(i,2) +vrho(2,j)
        exc(i) = exc(i) +exc0(j) *( rho(1,j) +rho(2,j) )

        dFc_daa(i) = v2rho2(1,j)
        dFc_dab(i) = v2rho2(2,j)
        dFc_dbb(i) = v2rho2(3,j)

#ifdef KXC_LXC_AVAILABLE
        dFc_daaa(i) = v3rho3(1,j)
        dFc_daab(i) = v3rho3(2,j)
        dFc_dabb(i) = v3rho3(3,j)
        dFc_dbbb(i) = v3rho3(4,j)
#endif
     End Do
  End if

  deallocate( rho );
  deallocate( vrho )
  deallocate( v2rho2 )
#ifdef KXC_LXC_AVAILABLE
  deallocate( v3rho3 )
  deallocate( v4rho4 )
#endif
  deallocate( exc0 )

end subroutine cr_lda_paw_libxc

subroutine set_hyrid_params_exch
  use m_Const_Parameters,  only : DP, ON
  use m_Control_Parameters,  only : xc_info_exch, xc_name_exch, xc_func_exch, &
       &                            alpha_exx, omega_exx, omega_exx_pbe, &
       &                            sw_screened_exchange
  use m_Parallelization,     only : mype

  use xc_f03_lib_m

  implicit none

  integer :: num, i
  real(kind=DP) :: c1
  real(kind=DP), allocatable :: params(:)

  alpha_exx = xc_f03_hyb_exx_coef( xc_func_exch )

  num = xc_f03_func_info_get_n_ext_params( xc_info_exch )
  if ( num <= 0 ) return

  allocate( params(num) )
  Do i=0, num -1
     c1 = xc_f03_func_info_get_ext_params_default_value( xc_info_exch, i )
     params(i+1) = c1
  End Do
  
  deallocate( params )

end subroutine set_hyrid_params_exch

subroutine set_hyrid_params_corr
  use m_Const_Parameters,  only : DP, ON
  use m_Control_Parameters,  only : xc_info_corr, xc_name_corr, xc_func_corr, &
       &                            alpha_exx, omega_exx, omega_exx_pbe, &
       &                            sw_screened_exchange
  use m_Files,              only : nfout

  use xc_f03_lib_m

  implicit none

  integer :: num, i
  real(kind=DP) :: c1
  real(kind=DP), allocatable :: params(:)

  alpha_exx = 0.25d0
  c1 = xc_f03_hyb_exx_coef( xc_func_corr )
  if ( c1 /= 0.0d0 ) alpha_exx = c1

  write(nfout,*) " Exx mixing parameter (trial)    = ", alpha_exx
  write(nfout,*)

  if ( xc_name_corr == "hyb_gga_xc_hse03" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse06" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse12" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse_sol"  &
       &         .or. xc_name_corr == "hyb_gga_xc_hjs_pbe"   &
       &         .or. xc_name_corr == "hyb_gga_xc_hjs_pbe_sol"   ) then
     sw_screened_exchange = ON
  endif

  num = xc_f03_func_info_get_n_ext_params( xc_info_corr )
  if ( num <= 0 ) return

  allocate( params(num) )
  Do i=0, num -1
     c1 = xc_f03_func_info_get_ext_params_default_value( xc_info_corr, i )
     params(i+1) = c1
  End Do
  
  if ( xc_name_corr == "hyb_gga_xc_hse03" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse06" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse12" &
       &         .or. xc_name_corr == "hyb_gga_xc_hse_sol"  &
       &         .or. xc_name_corr == "hyb_gga_xc_hjs_pbe"   &
       &         .or. xc_name_corr == "hyb_gga_xc_hjs_pbe_sol"   ) then
     alpha_exx     = params(1)        ! mixing parameter
     omega_exx     = params(2)        ! screening parameter for HF
     omega_exx_pbe = params(3)        ! screening parameter for PBE
!
     write(nfout,*) " Exx mixing parameter      = ", params(1)
     write(nfout,*) " screening parameter (HF)  = ", params(2)
     write(nfout,*) " screening parameter (DFT) = ", params(3)
     write(nfout,*)
  endif

  deallocate( params )

end subroutine set_hyrid_params_corr
