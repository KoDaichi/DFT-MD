!================================================
!  Software name : STM
!  Subroutine(s) : rd_WFs_doFFT_and_solve_eq, alloc_afft_etc, dealloc_afft_etc,
!                  print_ck_etc,print_VLC_Rspace, sumup_afft_to, print_dz_etc,
!                  print_summative_charges
!  Author(s)     : Koichi Kato and Takahiro Yamasaki (June 7, 2004)
!
!  FURTHER MODIFICATION: Junichiro  Koga (June 24, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

subroutine rd_WFs_doFFT_and_solve_eq
! $Id: rd_WFs_doFFT_and_solve_eq.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_Const_Parameters,     only : DP, BUCS, CARTS &
       &                           , HIGHER, LOWER, IN_BETWEEN, INVERSE, DIRECT &
       &                           , Hartree, NEGATIVE
  use m_Control_Parameters,   only : e1, e2, izi, izf,rini,rfin,nfin,erlmt &
       &                           , n_fc,n_fc_i
  use m_Charge_File,          only : natm2,natm,x_origin,y_origin,z_origin &
       &                           , fft_param,cell1,cell2,cell3 &
       &                           , iatomtype,ival,atom_pos
  use m_Kpoints, only :              kv3, nspin, kweight
  use m_Electronic_Structure, only : neordr,eko,efermi &
       &                           , m_ES_rd_WFs, m_ES_rd_VLCs &
       &                           , m_ES_WF_in_Rspace_fine_mesh &
       &                           , m_ES_VLC_in_Rspace_fine_mesh
  use m_FFT,                  only : nfft, nfftp, nfftpf, nlpf, fft_box_size_fine &
       &                           , m_FFT_WF_2dFFT_fine, m_FFT_Vlocal_W_fine
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Crystal_Structure,    only : neg, zl
  use m_Files,                only : nfzaj, nfinp, nfout, nfchgu, nfchgd, nfvlc, nfchgu_p &
       &                           , nfchgd_p, nfvlcr
  use m_Kpoints,              only : vkxyz
  use m_PlaneWaveBasisSet,    only : m_pwBS_kinetic_energies
  use m_ArraySize_Parameters, only : kspin, kimg

  implicit none
  integer :: i, kg2d, nls, it, ig, nlfid, n,n0, nmp, nnp
  real(kind=DP) :: ertmp, ermx, rmix, dz, ck, rdz, dz2, errsum
  real(kind=DP), allocatable, dimension(:,:,:,:) :: sfft, tfft
  real(kind=DP), allocatable, dimension(:,:,:) :: c, y
  real(kind=DP), allocatable, dimension(:,:)   :: am, w
  real(kind=DP), allocatable, dimension(:) :: ekin
  real(kind=DP), allocatable, dimension(:) :: bfft, cfft
  real(kind=DP), allocatable, dimension(:) :: ertmp_t, domin
  real(kind=DP), allocatable, dimension(:)  :: afft
  integer                                   :: id_sname = -1
  integer                                   :: ik, is, ib, ibto, ikp

  nmp = fft_box_size_fine(2)
  nnp = fft_box_size_fine(3)

  dz = zl / (nlpf*(3-kimg))   ! unit length along z axis
  dz2 = dz*dz
  !rdz = 1./dz2            !  = ar(i), al(i), i=1,nls
  rdz = 0.5d0/dz2            !  = ar(i), al(i), i=1,nls
!!$  nls = nlpf - izi + 2    ! mesh number for WF tail
  nls = izf - izi + 2    ! mesh number for WF tail
  nlfid = fft_box_size_fine(0)
  kg2d = nmp * nnp        ! 2D Gpoints
  e1 = e1 / Hartree       !eV -> Hartree unit
  e2 = e2 / Hartree       !eV -> Hartree unit

  call alloc_afft_etc  !-(contained here)
  call print_dz_etc       ! -(contained here)
!
  call m_ES_rd_VLCs(nfvlc)               !-(m_E.S.)
  sfft = 0.d0 ; tfft = 0.d0
  Kpoints: do ik = 1, kv3-nspin+1, nspin

     do is = 1, nspin
        print '(" spin = ", i3)', is

        call m_ES_VLC_in_Rspace_fine_mesh(is,bfft) !-(m_E.S.) \Vlc(G) -> \Vlc(r)
        call print_VLC_Rspace(is)             !-(contained here)
     
!!$        cycle
!        ikp = is + kspin*(ik-1)
        ikp = ik+(is-1)
        print '(" ik = ", i3, " ikp = ",i3)', ik, ikp
        call m_pwBS_kinetic_energies(ikp,vkxyz,ekin,nmp,nnp) !-(m_P.W.B.S)
        call m_ES_rd_WFs(nfzaj)
        write(6, *) 'efermi',efermi
        do ib = 1, neg
           ibto = neordr(ib,ikp)
           print '(" ** ib, ibto = ",2i6)',ib, ibto
           if(eko(ibto,ikp) > efermi + e2 ) exit
           if(eko(ibto,ikp) < efermi + e1 ) cycle

           ck = vkxyz(ikp,n_fc(2),CARTS)**2+vkxyz(ikp,n_fc(3),CARTS)**2 &
!!$           ck = vkxyz(ikp,n_fc(2),BUCS)**2+vkxyz(ikp,n_fc(3),BUCS)**2 &
                & - eko(ibto,ikp)
           if(ck > 0.) then
              ck = sqrt(ck)
           else
              ck = NEGATIVE
           end if
           call print_ck_etc
!!$#ifdef _DEBUG_WRITE_
!!$           call print_ck_etc
!!$#endif

           call m_ES_WF_in_Rspace_fine_mesh(ikp,ibto,afft)  !-(m_E.S.) \Psi(G) -> \Psi(r)
           call sumup_afft_to(tfft(:,:,:,is), ikp)

#ifdef _DEBUG_WRITE_
           print '(" m_ES_WF_in_Rspace_fine_mesh, ib = ",i3)',ib
           print *,'Wave Functions afft(first row in R-space)'
           do i = 1, fft_box_size_fine(1)/2 + 1
              n = 2*i - 1
              write(nfout,9001) afft(n), afft(n+1)
           end do
#endif
!
! convolution of local potential and WF start
!
           it = 1
1000       continue

           call m_FFT_Vlocal_W_fine (afft,bfft,cfft) !-(m_FFT)  cfft<-VLC(r)*Psi(r)
           call m_FFT_WF_2dFFT_fine(cfft,DIRECT)     !-(m_FFT) \Vlc*WF(r) -> \Vlc*WF(G(2D)&r(1D))

           call m_FFT_WF_2dFFT_fine(afft,DIRECT)     !-(m_FFT) \WF(r) -> \WF(G(2D)&r(1D))

!!$#ifdef _DEBUG_WRITE_
           if( it .lt. 3 ) then
              write(nfout,*) 'Convolution of WF & Potentials FFT-Backward cfft(first row)'
              print *,'no1 ekin =',ekin(1)
              do i = izi, nlpf
                 n = 2*i - 1
                 write(nfout,'(2d20.8)') cfft(n), cfft(n+1)
              end do
           end if
!!$#endif

           call tstatc0_begin('rd_WFs_doFFT_and_solve_eq_core ',id_sname)
!    print *,'finished \WF(r) -> \WF(G(2D)&r(1D))'
!
! diagonalization for tail WF start
!
           rmix = what_rmix(rini,rfin,nfin,it)

           am = 0. ; c = 0.; y = 0.; w = 0.

           do i = 2, nls
!!$              ar(i) = 1.0/dz**2 ;  al(i) = 1.0/dz**2
              do ig = 1, kg2d
                 n  = 1 + 2 * (izi+i-1-1) + kimg*nlfid*(ig-1)
                 !                     nlfid = fft_box_size_fine(0)
!                 am(ig,i) = eko(ibto, ikp) - ekin(ig) - 2./dz**2 
                 am(ig,i) = eko(ibto, ikp) - ekin(ig) - 1./dz**2 
                 c(ig,i,1) = cfft(n)
                 c(ig,i,2) = cfft(n+1)
              end do
           end do

           do ig = 1, kg2d
              n  = 1 + 2*(izi+2-1-1) + kimg*nlfid*(ig-1)
              n0 = 1 + 2*(izi+1-1-1) + kimg*nlfid*(ig-1)
              c(ig,2, 1) = cfft(n)   - afft(n0)*rdz
              c(ig,2, 2) = cfft(n+1) - afft(n0+1)*rdz
              c(ig,nls, 1) = 0.d0
              c(ig,nls, 2) = 0.d0
           end do

!!$           if(ck > 0) then
!!$              do ig = 1, kg2d
!!$                 am(ig,nls) = (dz*ck + 2.)/((dz*ck - 2.)*dz**2)
!!$              end do
!!$           end if

              ! -- forward substitution -->
           do ig = 1, kg2d
!!$              l(ig,2) = ar(2) / am(ig,2) 
              w(ig,2) = rdz / (am(ig,2))
              y(ig,2,1) = c(ig,2,1)/am(ig,2)
              y(ig,2,2) = c(ig,2,2)/am(ig,2)
           end do

           do i = 3, nls - 1
              do ig = 1, kg2d
                 w(ig,i) = rdz /  (am(ig,i) - rdz * w(ig,i-1))
                 y(ig,i,1) = (c(ig,i,1) - rdz*y(ig,i-1,1))*w(ig,i)/rdz
                 y(ig,i,2) = (c(ig,i,2) - rdz*y(ig,i-1,2))*w(ig,i)/rdz
              end do
           end do

           do ig = 1, kg2d
              y(ig,nls,1) = (c(ig,nls,1)-rdz*y(ig,nls-1,1)) &
                   &        /(am(ig,nls) - rdz*w(ig,nls-1))
              y(ig,nls,2) = (c(ig,nls,2)-rdz*y(ig,nls-1,2)) &
                   &        /(am(ig,nls) - rdz*w(ig,nls-1))
              ! <--

              ! -- backward substitution -->
!!$              c(ig,nls, 1) = y(ig,nls,1)/(d(ig,nls)**2)
!!$              c(ig,nls, 2) = y(ig,nls,2)/(d(ig,nls)**2)
              c(ig,nls, 1) = y(ig,nls,1)
              c(ig,nls, 2) = y(ig,nls,2)
           end do

           do i = nls-1, 2, -1
              do ig = 1, kg2d
                 c(ig,i,1) = y(ig,i,1) - w(ig,i) * c(ig,i+1,1)
                 c(ig,i,2) = y(ig,i,2) - w(ig,i) * c(ig,i+1,2)
              end do
           end do
              ! <--

           if(it.eq.1) then
              print '( " c (1:nls) (R|I)")'
              do i = 1, nls
                 print '(2d12.4)', c(1,i,1),c(1,i,2)
              end do
           end if

           ermx = 0.
           errsum = 0.
           do ig = 1, kg2d
             do i = 2, nls
                n = 1 + 2*(izi+i-1-1) + kimg*nlfid*(ig-1)
                if(sqrt(afft(n)**2 + afft(n+1)**2) > 1.e-8) then
!!$                if(sqrt(afft(n)**2 + afft(n+1)**2) < 1.e-8) cycle
                   ertmp =  ((c(ig,i,1)-afft(n))**2 &
                        &  + (c(ig,i,2)-afft(n+1))**2) &
                        &   / (afft(n)**2 + afft(n+1)**2)
                   errsum = errsum + ertmp
                   if ( ertmp > ermx) ermx = ertmp
                   afft(n)   = (1-rmix)*afft(n)   + rmix*c(ig,i,1)
                   afft(n+1) = (1-rmix)*afft(n+1) + rmix*c(ig,i,2)
                end if
              end do
           end do
           ermx = dsqrt(ermx)

           print '(" ( ",i4," ) rmix = ",f6.2," ||d WF|| (errsum) = " &
                & ,d12.4,"  (",d12.4,")")' ,it, rmix, ermx,errsum
           it = it + 1

!       if( ermx .le. erlmt ) then
!         write(nfout,*) 'Wave Functions afft(1st row in 2g-space) recalculated'
!         print *,'no1 ekin =',ekin(1)
!         write(nfout,*) (afft(i), i=1,nlpfd+2)
!       end if

           call tstatc0_end(id_sname)

           call m_FFT_WF_2dFFT_fine(afft,INVERSE)  !-(m_FFT)  \WF(G(2D)&r(1D)) -> \WF(r)

           if(it > 10*nfin) stop ' !! it > 10*nfin, but no-conversion'
           if(ermx > erlmt ) go to 1000

           print '(" *** Wave Function Converged!! ***,  it = ",i6)',it

#ifdef _DEBUG_WRITE_
           write(nfout,*) 'Wave Functions afft(first row in R-space) recalculated'
           write(nfout,*) (afft(i), i=1,2*(izf+1))
           write(nfout,*) 'Electron Densities afft(first row) recalculated'
           write(nfout,*) (afft(2*i-1)**2 + afft(2*i)**2, i=1, izf + 1)
#endif

           call sumup_afft_to(sfft(:,:,:,is), ikp)
        end do
     end do
  end do Kpoints
  do is=1,nspin
     call print_summative_charges(nfchgu,nfchgd,sfft(:,:,:,is),e1,e2)
     call print_summative_charges(nfchgu_p,nfchgd_p,tfft(:,:,:,is),e1,e2)
  enddo
  call dealloc_afft_etc
contains
!!$  subroutine solve_txdsml_equation(ibto,ikp,kg2d,ck)
!!$    integer,       intent(in) :: ibto, ikp, kg2d
!!$    real(kind=DP), intent(in) :: ck
!!$
!!$    integer :: i, ig
!!$
!!$    do ig = 1, kg2d
!!$       w(ig, 2) = ar(2) / am(ig, 2) 
!!$       g(ig,2,1) = c(ig,2,1)/am(ig,2)
!!$       g(ig,2,2) = c(ig,2,2)/am(ig,2)
!!$    end do
!!$
!!$    do i = 3, nls - 1
!!$       do ig = 1, kg2d
!!$          w(ig,i) = ar(i) / (am(ig,i) - al(i) * w(ig,i-1))
!!$          g(ig,i,1) = (c(ig,i,1) - al(i)*g(ig,i-1,1))*w(ig,i)/ar(i)
!!$          g(ig,i,2) = (c(ig,i,2) - al(i)*g(ig,i-1,2))*w(ig,i)/ar(i)
!!$       end do
!!$    end do
!!$    do ig = 1, kg2d
!!$       g(ig,nls,1) = (c(ig,nls,1) - al(nls)*g(ig,nls-1,1)) &
!!$            & /(am(ig,nls) - al(nls)*w(ig,nls-1))
!!$       g(ig,nls,2) = (c(ig,nls,2) - al(nls)*g(ig,nls-1,2)) &
!!$            & /(am(ig,nls) - al(nls)*w(ig,nls-1))
!!$
!!$       c(ig,nls,1) = g(ig,nls,1)
!!$       c(ig,nls,2) = g(ig,nls,2)
!!$    end do
!!$    do i = nls-1, 2, -1
!!$       do ig = 1, kg2d
!!$          c(ig,i,1) = g(ig,i,1) - w(ig,i) * c(ig,i+1,1)
!!$          c(ig,i,2) = g(ig,i,2) - w(ig,i) * c(ig,i+1,2)
!!$       end do
!!$    end do
!!$
!!$  end subroutine solve_txdsml_equation

  subroutine alloc_afft_etc
    allocate(afft(nfftpf))
    allocate(bfft(nfftpf))
    allocate(cfft(nfftpf))
    allocate(sfft(nmp,nnp,nlpf+1,nspin))
    allocate(tfft(nmp,nnp,nlpf+1,nspin))
    allocate(c(nmp*nnp+1,nls, 2))
    allocate(y(nmp*nnp+1,nls, 2))
!!$    allocate(c(izf+2, 2))
!!$    allocate(g(izf+2, 2))
    allocate(ekin(nmp*nnp+1))
!!$    allocate(al(nls))
    allocate(am(nmp*nnp+1, nls))
    allocate(w(nmp*nnp+1,nls))
!!$    allocate(ar(nls))
!!$    allocate(l(nmp*nnp+1, nls))
    allocate(ertmp_t(nmp*nnp))
    allocate(domin(nmp*nnp))
  end subroutine alloc_afft_etc

  subroutine dealloc_afft_etc
    deallocate(afft)
    deallocate(bfft)
    deallocate(cfft)
    deallocate(sfft)
    deallocate(tfft)
    deallocate(c)
    deallocate(y)
    deallocate(ekin)
!!$    deallocate(al)
    deallocate(am)
    deallocate(w)
!!$    deallocate(ar)
!!$    deallocate(l)
    deallocate(ertmp_t)
    deallocate(domin)
  end subroutine dealloc_afft_etc
    
  subroutine print_ck_etc
    print *,'ck =', ck
    print '(" ib = ",i3," ibto = ", i3)', ib, ibto
    print '(" eko(ibto,is + kspin*(ik-1)) = ",f10.6)',eko(ibto,ikp)
  end subroutine print_ck_etc
    
  subroutine print_VLC_Rspace(is)
    integer, intent(in) :: is
    integer :: i, n, idp,nlp,nmp,nnp,nlphf, j,k, ip(4),knew,inew,jnew,nlph, ixh,iyh,ipx(4),ipy(4)
    real(kind=DP), allocatable, dimension(:,:,:) :: wk

    print *,'finished m_ES_VLC_in_Rspace_fine_mesh'
    ixh = fft_box_size_fine(2)/2
    iyh = fft_box_size_fine(3)/2
    ipx(1) = 0;     ipy(1) = 0
    ipx(2) = ixh-1; ipy(2) = 0
    ipx(3) = 0;     ipy(3) = iyh-1
    ipx(4) = ixh-1; ipy(4) = iyh-1

    write(nfout,'(" kimg = ",i6," fft_box_size_fine(1) = ",i8)') kimg, fft_box_size_fine(1)
!!$    write(nfout,*) 'Local Potentials (first row in R-space)'
    write(nfout,*) 'Local Potentials ( 4 rows in R-space)'
    write(nfout,*) ' first coloumn:  ',ipx(1),ipy(1)
    write(nfout,*) ' second coloumn: ',ipx(2),ipy(2)
    write(nfout,*) ' third coloumn:  ',ipx(3),ipy(3)
    write(nfout,*) ' fourth coloumn: ',ipx(4),ipy(4)

!!$    do it = 1, 2
!!$       if(it == 1) then
!!$       else if(it == 2) then
!!$          write(nfout,*) 'Local Potentials (second row in R-space)'
!!$       end if
!!$    end do

    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1); nmp = fft_box_size_fine(2); nnp = fft_box_size_fine(3)
    if(kimg.eq.1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if
    
    do i = 1, 4
       ip(i) = (nlphf*nmp*ipy(i) + nlphf*ipx(i))*2
    end do
    if(kimg.eq.2) then
!!$       if(it == 2) ip = fft_box_size_fine(0)*2
       do i = 1, fft_box_size_fine(1)
          n = 2*i-1
!!$          write(nfout,'(i8,2d20.8)') i,bfft(ip+n),bfft(ip+n+1)
          write(nfout,'(i8,4d15.7)') i,bfft(ip(1)+n),bfft(ip(2)+n),bfft(ip(3)+n),bfft(ip(4)+n)
       end do
    else if(kimg.eq.1) then
!!$       if(it == 2) ip = fft_box_size_fine(0)
       do i = 1, fft_box_size_fine(1)/2
!!$          write(nfout,'(i8,2d20.8)') i,bfft(ip+2*i-1),bfft(ip+2*i)
          write(nfout,'(i8,4d15.7)') i,bfft(ip(1)+i),bfft(ip(2)+i),bfft(ip(3)+i),bfft(ip(4)+i)
       end do
    end if

    write(nfvlcr,*) ' spin = ', is
    allocate(wk(nlp,nmp,nnp))
    nlph = nlp
    if(kimg.eq.1) nlph = nlp/2
    do i = 1, nmp
       do j = 1, nnp
          do k = 1, nlp
             if(kimg.eq.1.and.k.gt.nlphf) then
                knew = idp -k
                jnew = nnp+2 -j
                inew = nmp+2 -i
                if(jnew.gt.nnp) then
                   jnew = jnew - nnp
                end if
                if(inew.gt.nmp) then
                   inew = inew - nmp
                end if
             else
                knew = k
                jnew = j
                inew = i
             end if
             ip(1) = nlphf*nmp*(jnew-1) + nlphf*(inew-1) + knew
             wk(k,i,j) = bfft(ip(1)*2-1)
          end do
       end do
    end do
    write(nfout,*) ' !D FFT cube mapping finished'
    write(nfvlcr,9001) nlp*nmp*nnp,nlp,nmp,nnp
    write(nfvlcr,9002) wk
 9001 format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')
 9002 format(6e12.4)
  end subroutine print_VLC_Rspace

  subroutine sumup_afft_to(xfft, ik)
!    real(kind=DP), intent(out), dimension(nmp,nnp,nlpf+1) :: xfft
    real(kind=DP), intent(inout), dimension(nmp,nnp,nlpf+1) :: xfft
    integer, intent(in) :: ik
    integer :: k, j, i, n
    do k = 1, nlpf
       do j = 1, nnp
          do i = 1, nmp
             n = 1 + 2*(k-1) + kimg*nlfid*((i-1) + nmp*(j-1))
             xfft(i,j,k) = xfft(i,j,k) + kweight(ik)*(afft(n)**2 + afft(n+1)**2)
          end do
       end do
    end do

!!$    do k = 1, nnp
!!$       do j = 1, nmp
!!$          do i = 1, nlpf + 1
!!$             n = 1 + 2*(i-1) + kimg*nlfid*((j-1) + nmp*(k-1))
!!$             sfft(i, j, k) = sfft(i, j, k) + afft(n)**2 + afft(n+1)**2
!!$          end do
!!$       end do
!!$    end do
  end subroutine sumup_afft_to

  real(kind=DP) function what_rmix(rini,rfin,nfin,it)
    real(kind=DP), intent(in) :: rini,rfin
    integer, intent(in)       :: nfin,it
    if( it <= nfin ) then
       what_rmix = rini + (rfin - rini) * (it-1)/nfin
    else
       what_rmix = rfin
    end if
!!$    print '(" Iteration Loop No, rmix = ",i8,f10.4)', it, rmix
!!$    write(nfout,*) 'Iteration Loop No, rmix = ', it, rmix
  end function what_rmix

  subroutine print_dz_etc
    print '(" -- dz -- ", f10.6)', dz
    print '(" -- nls, nlpf, izi, izf -- ", 4i6)', nls,nlpf,izi,izf
    print '(" -- rini, rfin, nfin -- ",2f10.6,i5)',rini,rfin,nfin
    print '(" -- erlmt -- ",f10.6)',erlmt
    print '(" unit length along z axis = ",f10.6)',dz
    print '(" Bias Voltage (V) = ", f10.6)', e1
  end subroutine print_dz_etc
    
  subroutine print_summative_charges(nfchgu,nfchgd,xfft,e1,e2)
    integer, intent(in) :: nfchgu, nfchgd
    integer :: i,j,k,scale
    integer, dimension(3) :: index_tmp
    real(kind=DP), intent(in) :: xfft(nmp,nnp,nlpf+1),e1,e2
    real(kind=DP), allocatable, dimension(:,:,:) :: xfft_wk
    integer :: nfchg
!
! Print out charges sumed at all k points
!
    if ( is .eq. 1) then
       nfchg = nfchgu
       print *,' -- nfchgu --'
    else if( is .eq. 2) then
       nfchg = nfchgd
       print *,' -- nfchgd --'
    end if


!      write(nfchg,9001) nlpf*nmp*nnp,nmp,nnp,nlpf, e1*Hartree, e2*Hartree
!      write(nfchg,501) xfft
!!$    do k = 1, nnp
!!$       do j = 1, nmp
!!$          do i = 1, nlpf + 1
!!$             write(nfchr,500) i, j, k, sfft(i,j,k)
!!$          end do
!!$       end do
!!$    end do

      write(nfchg,9001) nlpf*nmp*nnp,nmp,nnp,nlpf, e1*Hartree, e2*Hartree
      if ( is .eq. 1) then
         write(nfchg,*) ' -- nfchgu --'
      else if( is .eq. 2) then
         write(nfchg,*) ' -- nfchgd --'
      end if

      if ( kimg .eq. 1 ) then
        scale = 2
      else
        scale = 1
      endif

      write(nfchg,'(i6,3f12.6)') natm, x_origin,y_origin,z_origin

      index_tmp(n_fc(1)) = nlpf
      index_tmp(n_fc(2)) = nmp
      index_tmp(n_fc(3)) = nnp

      write(nfchg,'(i6,3f12.6)') index_tmp(1),cell1(1:3)
      write(nfchg,'(i6,3f12.6)') index_tmp(2),cell2(1:3)
      write(nfchg,'(i6,3f12.6)') index_tmp(3),cell3(1:3)

      do i=1,natm
        write(nfchg,'(i6,4f12.6)') iatomtype(i),ival(i),atom_pos(i,1:3)
      enddo
        
      allocate(xfft_wk(index_tmp(3),index_tmp(2),index_tmp(1)))

      do j=1,nnp
        do k=1,nmp
          do i=1,nlpf
            index_tmp(n_fc_i(1)) = j
	    index_tmp(n_fc_i(2)) = k
	    index_tmp(n_fc_i(3)) = i
	    xfft_wk(index_tmp(1),index_tmp(2),index_tmp(3)) = xfft(k,j,i)
	  enddo
	enddo
      enddo

      write(nfchg,501) xfft_wk
      deallocate(xfft_wk)
     
    print *,'finished write nfchg'

 9001 format(' CHARGE DENSITY NE = ',i8,'(',3i5,' )  Energy-Range ',f8.4,' : ',f8.4,' (eV)')
500   format(3I4,1pe12.4)
501   format(6e16.8)

  end subroutine print_summative_charges

end subroutine rd_WFs_doFFT_and_solve_eq
