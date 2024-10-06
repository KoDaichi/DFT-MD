!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.00)
!
!  MODULE: m_FFT
!
!  AUTHOR(S): T. Yamasaki, K. Betsuyaku,   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
!                                   , March/14/2005, April/10/2007
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
!
  integer :: idp, idp2, idp3, nlp, nmp, nnp
  integer, private      :: istat = 0
  integer :: sw_avoiding_odd_fftbox = ON
  integer :: sw_zero_padding = ON
  real(kind=DP),private,pointer,dimension(:)               :: cfft
  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1,cw2,cw3
  real(kind=DP),private,allocatable,dimension(:)           :: trigs_WF, trigs_pWF, trigs_CD
  integer,      private,allocatable,dimension(:)           :: ifax_WF, ifax_pWF, ifax_CD
  real(kind=DP),private,allocatable, dimension(:)          :: ftw
  ! ------- Positron start
  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1_pstrn,cw2_pstrn,cw3_pstrn
  ! ------- Positron end

!  include 'mpif.h'
contains 
  subroutine fft_WFCD_work_alloc
    integer :: nnl,nnm,nnn
! ----> work arrays of fft for wave functions
    nnl  = fft_box_size_WF(1,1)
    nnm  = fft_box_size_WF(2,1);  nnn  = fft_box_size_WF(3,1)
    if(ipri >= 1) write(nfout,'(" !!allocation of cw1,cw2,cw3 (fft_WFCD_work_alloc)")')
    if(kimg == 1) then
       allocate(trigs_WF(nnl+2*(nnm+nnn)), stat=istat)
       allocate(ifax_WF(60), stat=istat)
       if(ipri >= 1) write(nfout,'(" !! WF_ASLFFT <<m_FFT.fft_WFCD_work_alloc>>")')
    else if(kimg == 2) then
       allocate(trigs_WF(2*(nnl+nnm+nnn)), stat=istat)
       allocate(ifax_WF(60), stat=istat)
       if(ipri >= 1) write(nfout,'(" !! WF_ASLFFT")')
    end if
! ------- Positron start
! ----> and work arrays of fft for positron wave functions
    if(sw_positron /= OFF) then
       nnl  = fft_box_size_pWF(1,1)
       nnm  = fft_box_size_pWF(2,1);  nnn  = fft_box_size_pWF(3,1)
       if(ipri >= 1) write(nfout,'(" !!allocation of cw1_pstrn,cw2_pstrn,cw3_pstrn (fft_WFCD_work_alloc)")')
       if(kimg == 1) then
          allocate(trigs_pWF(nnl+2*(nnm+nnn)), stat=istat)
          allocate(ifax_pWF(60), stat=istat)
       else if(kimg == 2) then
          allocate(trigs_pWF(2*(nnl+nnm+nnn)), stat=istat)
          allocate(ifax_pWF(60), stat=istat)
       end if
    end if
! ------- Positron end

  end subroutine fft_WFCD_work_alloc

  subroutine m_FFT_alloc_WF_work
  end subroutine m_FFT_alloc_WF_work

  subroutine m_FFT_alloc_pWF_work()
  end subroutine m_FFT_alloc_pWF_work

  subroutine m_FFT_dealloc_WF_work
  end subroutine m_FFT_dealloc_WF_work

  subroutine m_FFT_alloc_CD_box
  end subroutine m_FFT_alloc_CD_box

  subroutine m_FFT_dealloc_CD_box
  end subroutine m_FFT_dealloc_CD_box

  subroutine m_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset
!!$    integer :: nnl,nnm,nnn
    integer :: nx,ny,nz, ic, isize

    integer :: nfft_t
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup ',id_sname)


! --- fft_box_size_WF, fft_box_size_pWF ---
    if(ipri >= 1) write(nfout,*) '!FFT_WF = ASLFFT <<m_FFT_setup>>'
    if(inversion_symmetry == ON) then  ! kimg == 1
       if(mod(fft_box_size_WF(1,1),4) == 2) then
          fft_box_size_WF(1,0) = fft_box_size_WF(1,1) + 4
       else
          fft_box_size_WF(1,0) = fft_box_size_WF(1,1) + 2
       end if
       fft_box_size_WF(2:3,0) = fft_box_size_WF(2:3,1) + 1
       if(sw_positron /= OFF) then
          if(mod(fft_box_size_pWF(1,1),4) == 2) then
             fft_box_size_pWF(1,0)   = fft_box_size_pWF(1,1) + 4
          else
             fft_box_size_pWF(1,0)   = fft_box_size_pWF(1,1) + 2
          end if
          fft_box_size_pWF(2:3,0) = fft_box_size_pWF(2:3,1) + 1
       end if
    else if(inversion_symmetry == OFF) then ! kimg == 2
       fft_box_size_WF(1:3,0) = fft_box_size_WF(1:3,1) + 1
       if(sw_positron /= OFF) then
          fft_box_size_pWF(1:3,0) = fft_box_size_pWF(1:3,1) + 1
       end if
    endif
! --- fft_box_size_CD ---
    if(ipri >= 1) write(nfout,*) '!FFT_CD = MPIFFT <<m_FFT_setup>>'
    if(inversion_symmetry == ON) then  ! kimg == 1
       if(mod(fft_box_size_CD(1,1),4) == 2) then
          fft_box_size_CD_nonpara(1,0) = fft_box_size_CD(1,1) + 4
       else
          fft_box_size_CD_nonpara(1,0) = fft_box_size_CD(1,1) + 2
       end if
       fft_box_size_CD_nonpara(2:3,0) = fft_box_size_CD(2:3,1) + 1
    else if(inversion_symmetry == OFF) then ! kimg == 2
       fft_box_size_CD_nonpara(1:3,0) = fft_box_size_CD(1:3,1) + 1
    endif
!!$    if(inversion_symmetry == ON) then  ! kimg == 1
!!$       ipad = 2
!!$    else if(inversion_symmetry == OFF) then ! kimg == 2
!!$       ipad = 1
!!$    end if
!!$    fft_box_size_CD_nonpara(1,0)   = fft_box_size_CD(1,1) + ipad
!!$    fft_box_size_CD_nonpara(2:3,0) = fft_box_size_CD(2:3,1)
    nfftp_nonpara  = product(fft_box_size_CD_nonpara(1:3,0)) * (2-inversion_symmetry)

    nx = fft_box_size_CD(1,1); ny = fft_box_size_CD(2,1); nz = fft_box_size_CD(3,1)
    do ic = 1, 2
       if(ic == 1) isize = npes
       if(ic == 2) isize = npes_cdfft
       ny_d =  ny/isize
       nz_d =  nz/isize
       if(mod( ny,isize) .ne. 0) ny_d = ny_d + 1
       if(mod( nz,isize) .ne. 0) nz_d = nz_d + 1
       ly_d = ny_d;   lz_d = nz_d
       if(inversion_symmetry == ON) then  ! kimg == 1
          lx = nx + 2
       else if(inversion_symmetry == OFF) then ! kimg == 2
          lx = nx
          if(mod(lx,2) == 0)   lx = nx + ldx
       end if

       if(mod(ly_d,2) == 0) ly_d = ny_d + ldy
       if(mod(lz_d,2) == 0) lz_d = nz_d + ldz

       ly   = ly_d*isize
       lz   = lz_d*isize
       if(ic == 1) then
          fft_box_size_CD(1,0) = lx
          fft_box_size_CD(2,0) = ly
          fft_box_size_CD(3,0) = lz
          ny_ds = ny_d;          nz_ds = nz_d;          ly_ds = ly_d
          lys   = ly;            lz_ds = lz_d;          lzs   = lz
          lxs   = lx
       else
          fft_box_size_CD_c(1,0) = lx
          fft_box_size_CD_c(2,0) = ly
          fft_box_size_CD_c(3,0) = lz
       end if
    end do
    if(ipri >= 1) then
       write(nfout,'(" nx,  ny,  nz,  ny_d, nz_d  = ",5i8)') nx,ny,nz,ny_d,nz_d
       write(nfout,'("      ly,  lz,  ly_d, lz_d  = ",8x,4i8)') ly,lz,ly_d,lz_d
       write(nfout,'("                ny_ds,nz_ds = ",24x,2i8)') ny_ds, nz_ds
       write(nfout,'("      lys, lzs, ly_ds,lz_ds = ",8x,4i8)') lys,lzs,ly_ds,lz_ds
    end if

    nfft =   product(fft_box_size_WF(1:3,0)) * (2-inversion_symmetry)
    nfftp  = product(fft_box_size_CD_c(1:3,0)) * (2-inversion_symmetry)
    nfftps = product(fft_box_size_CD(1:3,0)) * (2-inversion_symmetry)

    if(ipri >= 1) call wd_FFTboxsizes(nfout)
    if(sw_positron /= OFF) &
         & nfft_pstrn = product(fft_box_size_pWF(1:3,0))*(2-inversion_symmetry)
    nfft_t = nfft
    if(sw_positron /= OFF .and. nfft_t < nfft_pstrn) nfft_t = nfft_pstrn

    if(.not. paramset) then
       call fft_WFCD_work_alloc     ! <cw[123],cw[123]_pstrn,wlp,wmp,wnp> are allocated

! Initialization of the Wave-Function FFT
       allocate(cfft(nfft_t), stat=istat) ! cfft is used only for initiallization
       if(istat /= 0) then
          if(ipri >= 1) then
             write(nfout,*) 'Allocation error for cfft in sub. m_FFT_setup'
             write(nfout,*) 'stat =', istat, 'nfft =', nfft
          end if
          stop
       end if

       call init_fft_coefficients_arrays_WF()

       deallocate(cfft,stat=istat)
       if(istat /= 0 ) then
          if(ipri >= 1) then
             write(nfout,*) 'Deallocation error for cfft in sub. m_FFT_setup'
             write(nfout,*) 'stat =', istat
          end if
          stop
       end if
    endif

    call tstatc0_end(id_sname)

  contains
    subroutine init_fft_coefficients_arrays_WF()
      use m_Parallelization,   only : itask
      integer :: id, nl, nm, nn, ierr
      integer :: id_p, nl_p, nm_p, nn_p
      real(kind=DP) :: ftw(nfft+nfft_pstrn)
      integer :: id2, id3, id_p2, id_p3

      id = fft_box_size_WF(1,0)
      id2 = fft_box_size_WF(2,0)
      id3 = fft_box_size_WF(3,0)
      nl = fft_box_size_WF(1,1)
      nm = fft_box_size_WF(2,1)
      nn = fft_box_size_WF(3,1)
      if(sw_positron /= OFF) then
         id_p = fft_box_size_pWF(1,0)
         id_p2 = fft_box_size_pWF(2,0)
         id_p3 = fft_box_size_pWF(3,0)
         nl_p = fft_box_size_pWF(1,1)
         nm_p = fft_box_size_pWF(2,1)
         nn_p = fft_box_size_pWF(3,1)
      end if

      if(kimg == 1) then
         call qfr3fb(nl,nm,nn,cfft,id,id2,id3,0,ifax_WF,trigs_WF,ftw,itask,ierr)
         if(sw_positron /= OFF) then
            call qfr3fb(nl_p,nm_p,nn_p,cfft,id_p,id_p2,id_p3,0,ifax_pWF,trigs_pWF,ftw,itask,ierr)
         endif
      else
         call hfc3fb(nl,nm,nn,cfft,id,id2,id3,0,ifax_WF,trigs_WF,ftw,itask,ierr)
         if(sw_positron /= OFF) then
            call hfc3fb(nl_p,nm_p,nn_p,cfft,id_p,id_p2,id_p3,0,ifax_pWF,trigs_pWF,ftw,itask,ierr)
         end if
      endif
    end subroutine init_fft_coefficients_arrays_WF
  end subroutine m_FFT_setup

  subroutine CDFFT_setup()
     CD_setup_is_done = YES
  end subroutine CDFFT_setup

  subroutine m_FFT_WF(electron_or_positron,nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    use m_Parallelization,   only : itask
    integer, intent(in)          :: electron_or_positron
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch
    complex(kind=CMPLDP),pointer,dimension(:) :: cw1_t, cw2_t, cw3_t
    integer :: id, nl, nm, nn, ier1 = 0
    real(kind=DP) :: ftw(nfft)
    integer, dimension(2) :: flag_asl_rfft      = (/-1,  1/)
    integer, dimension(2) :: flag_asl_cfft      = (/ 1, -1/)
    integer :: id2, id3
    integer :: isw
    integer :: id_sname = -1

    if(electron_or_positron == ELECTRON) then
       id = fft_box_size_WF(1,0)
       id2 = fft_box_size_WF(2,0)
       id3 = fft_box_size_WF(3,0)
       nl = fft_box_size_WF(1,1)
       nm = fft_box_size_WF(2,1)
       nn = fft_box_size_WF(3,1)
       cw1_t=>cw1; cw2_t=>cw2; cw3_t=>cw3
    else if(electron_or_positron == POSITRON) then
       id = fft_box_size_pWF(1,0)
       id2 = fft_box_size_pWF(2,0)
       id3 = fft_box_size_pWF(3,0)
       nl = fft_box_size_pWF(1,1)
       nm = fft_box_size_pWF(2,1)
       nn = fft_box_size_pWF(3,1)
       cw1_t=>cw1_pstrn; cw2_t=>cw2_pstrn; cw3_t=>cw3_pstrn
    end if

    if(electron_or_positron == ELECTRON) then
       if(kimg == 1) then
          isw = -flag_asl_rfft(inverse_or_direct)
!!$          if(ipri >= 2) write(nfout,*) 'size(afft)=',size(afft)
!!$          if(ipri >= 2) write(nfout,*) 'isw=',isw
!!$          if(ipri >= 2) write(nfout,*) 'nl,nm,nn,id,id2,id3=',nl,nm,nn,id,id2,id3
          call qfr3bf(nl,nm,nn,afft,id,id2,id3,isw,ifax_WF,trigs_WF,ftw,itask,ier1)
       else if(kimg == 2) then
          isw = -flag_asl_cfft(inverse_or_direct)
          call hfc3bf(nl,nm,nn,afft,id,id2,id3,isw,ifax_WF,trigs_WF,ftw,itask,ier1)
       end if
    else if(electron_or_positron == POSITRON) then
       if(kimg == 1) then
          isw = -flag_asl_rfft(inverse_or_direct)
          call qfr3bf(nl,nm,nn,afft,id,id2,id3,isw,ifax_pWF,trigs_pWF,ftw,itask,ier1)
       else if(kimg == 2) then
          isw = flag_asl_cfft(inverse_or_direct)
          call hfc3bf(nl,nm,nn,afft,id,id2,id3,isw,ifax_pWF,trigs_pWF,ftw,itask,ier1)
       end if
    end if

    call sxfft_zero_padding(afft,nl,nm,nn,id,id2,id3)

    if (ier1 /= 0) then
       if(ipri>=1) write (nfout,*) ' !!(Error fft) ier1 = ',ier1
       stop
    end if
  end subroutine m_FFT_WF

 subroutine sxfft_zero_padding(afft,nl,nm,nn,nd1,nd2,nd3)
      use m_Const_Parameters, only : DP
      implicit none
      real(kind=DP), intent(out) :: afft(*)
      integer, intent(in) :: nl,nm,nn
      integer, intent(in) :: nd1,nd2,nd3
      integer :: i,j,k,ip
      integer :: nlh,ndh
      integer :: imax
      imax = nd1*nd2*nd3

      if(kimg==1) then
         nlh = nl/2
         ndh = nd1/2
         if(ipri >= 2) write(6,*) 'ASLFFT: nl,ndh=',nlh,ndh
         do j = nlh+2, ndh
            do i = 1, nd2*nn
               ip = ndh*(i-1) + j
               afft(ip*2-1) = 0.d0
               afft(ip*2) = 0.d0
            end do
         end do
         do j = nm+1, nd2
            do k = 1, nd3
               do i = 1, ndh
                  ip = i + ndh*(j-1) + ndh*nd2*(k-1)
                  afft(ip*2-1) = 0.d0
                  afft(ip*2) = 0.d0
               end do
            end do
         end do
         do k = nn+1, nd3
            do i = 1, ndh*nd2
               ip = i + ndh*nd2*(k-1)
               afft(ip*2-1) = 0.d0
               afft(ip*2) = 0.d0
            end do
         end do

      else if(kimg==2) then
         do j = nl+1, nd1
            do i = 1, nd2*nn
               ip = nd1*(i-1) + j
               afft(ip*2-1) = 0.d0
               afft(ip*2) = 0.d0
            end do
         end do
         do j = nm+1, nd2
            do k = 1, nn
               do i = 1, nl
                  ip = i + nd1*(j-1) + nd1*nd2*(k-1)
                  afft(ip*2-1) = 0.d0
                  afft(ip*2) = 0.d0
               end do
            end do
         end do
         do k = nn+1, nd3
            do i = 1, nd1*nd2
               ip = i + nd1*nd2*(k-1)
               afft(ip*2-1) = 0.d0
               afft(ip*2) = 0.d0
            end do
         end do
      end if
 end subroutine sxfft_zero_padding

 subroutine fft_CD_inverse_core(afft_CD)
      use m_Parallelization,   only : itask
!!$      real(kind=DP), intent(inout) :: afft_cd(nfftp)
      real(kind=DP), intent(inout) :: afft_cd(nfftp_nonpara)
      real(kind=DP) :: ftw(nfftp_nonpara)
      integer :: idp, idp2, idp3, nlp, nmp, nnp
      integer :: isw, ier1

      idp  = fft_box_size_CD_nonpara(1,0)
      idp2 = fft_box_size_CD_nonpara(2,0)
      idp3 = fft_box_size_CD_nonpara(3,0)
      nlp  = fft_box_size_CD(1,1)
      nmp  = fft_box_size_CD(2,1)
      nnp  = fft_box_size_CD(3,1)

      if(.not.allocated(ifax_CD)) allocate(ifax_CD(60))
      if(kimg == 1) then
         if(.not.allocated(trigs_CD)) allocate(trigs_CD(nlp+2*(nmp+nnp)))
      else if(kimg == 2) then
         if(.not.allocated(trigs_CD)) allocate(trigs_CD(2*(nlp+nmp+nnp)))
      end if

      if(kimg == 1) then
         isw =  1
         call qfr3bf(nlp,nmp,nnp,afft_CD,idp,idp2,idp3,isw,ifax_CD,trigs_CD,ftw,itask,ier1)
      else if(kimg == 2) then
         isw = -1
         call hfc3bf(nlp,nmp,nnp,afft_CD,idp,idp2,idp3,isw,ifax_CD,trigs_CD,ftw,itask,ier1)
      end if
      call sxfft_zero_padding(afft_CD,nlp,nmp,nnp,idp,idp2,idp3)
 end subroutine fft_CD_inverse_core

 subroutine fft_CD_direct_core(afft_CD)
      use m_Parallelization,   only : itask
!!$      real(kind=DP), intent(inout) :: afft_cd(nfftp)
      real(kind=DP), intent(inout) :: afft_cd(nfftp_nonpara)
      real(kind=DP) :: ftw(nfftp_nonpara)
      integer :: isw, ier1

      idp  = fft_box_size_CD_nonpara(1,0)
      idp2 = fft_box_size_CD_nonpara(2,0)
      idp3 = fft_box_size_CD_nonpara(3,0)
      nlp  = fft_box_size_CD(1,1)
      nmp  = fft_box_size_CD(2,1)
      nnp  = fft_box_size_CD(3,1)

      if(.not.allocated(ifax_CD)) allocate(ifax_CD(60))
      if(kimg == 1) then
         if(.not.allocated(trigs_CD)) allocate(trigs_CD(nlp+2*(nmp+nnp)))
      else if(kimg == 2) then
         if(.not.allocated(trigs_CD)) allocate(trigs_CD(2*(nlp+nmp+nnp)))
      end if

      if(kimg == 1) then
         isw = -1
         call qfr3bf(nlp,nmp,nnp,afft_CD,idp,idp2,idp3,isw,ifax_CD,trigs_CD,ftw,itask,ier1)
      else if(kimg == 2) then
         isw =  1
         call hfc3bf(nlp,nmp,nnp,afft_CD,idp,idp2,idp3,isw,ifax_CD,trigs_CD,ftw,itask,ier1)
      end if
      call sxfft_zero_padding(afft_CD,nlp,nmp,nnp,idp,idp2,idp3)
 end subroutine fft_CD_direct_core
