!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: dxml_fft, mulfac
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
      subroutine dxml_fft(iop,id,nl,nm,nn,afft)
c $Id: dxmlfft.f 570 2017-04-21 20:34:50Z yamasaki $
      implicit real*8(a-h,o-z)
      include '/usr/include/DXMLDEF.FOR'
      real*8 afft(*)
      record /DXML_d_fft_structure_3D/ fft_struct_cd,fft_struct_wd
      record /DXML_z_fft_structure_3D/ fft_struct_cz,fft_struct_wz
      save fft_struct_cd,fft_struct_wd
      save fft_struct_cz,fft_struct_wz
      save nfft_cd,nfft_cz,nfft_wd,nfft_wz
      save nbox_cd,nbox_cz,nbox_wd,nbox_wz
c
c     IOP=1..DIRECT, IOP=2..INVERSE (REAL TRANSFORM)
c     IOP=1..INVERSE, IOP=2..DIRECT (COMPLEX TRANSFORM)
c
      entry cd_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call dfft_apply_3d('r','c','f',afft,afft,id,nm,
     $                        fft_struct_cd,1,1,1)
        else
          call dfft_apply_3d('c','r','b',afft,afft,id,nm,
     $                        fft_struct_cd,1,1,1)
c         scale=nfft_cd
c         call dscal(nbox_cd,scale,afft,1)
          call mulfac(afft,nbox_cd,nfft_cd)
        endif
      return
c
      entry cz_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call zfft_apply_3d('c','c','b',afft,afft,id,nm,
     $                        fft_struct_cz,1,1,1)
c         scale=nfft_cz
c         call dscal(nbox_cz,scale,afft,1)
          call mulfac(afft,nbox_cz,nfft_cz)
        else
          call zfft_apply_3d('c','c','f',afft,afft,id,nm,
     $                        fft_struct_cz,1,1,1)
        endif
      return
c
      entry wd_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call dfft_apply_3d('r','c','f',afft,afft,id,nm,
     $                        fft_struct_wd,1,1,1)
        else
          call dfft_apply_3d('c','r','b',afft,afft,id,nm,
     $                        fft_struct_wd,1,1,1)
c         scale=nfft_wd
c         call dscal(nbox_wd,scale,afft,1)
          call mulfac(afft,nbox_wd,nfft_wd)
        endif
      return
c
      entry wz_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call zfft_apply_3d('c','c','b',afft,afft,id,nm,
     $                        fft_struct_wz,1,1,1)
c         scale=nfft_wz
c         call dscal(nbox_wz,scale,afft,1)
          call mulfac(afft,nbox_wz,nfft_wz)
        else
          call zfft_apply_3d('c','c','f',afft,afft,id,nm,
     $                        fft_struct_wz,1,1,1)
        endif
      return
          
      entry cd_init(id,nl,nm,nn)
        call dfft_init_3d(nl,nm,nn,fft_struct_cd,.true.)
        nbox_cd=id*nm*nn
        nfft_cd=nl*nm*nn
        return
      entry cz_init(id,nl,nm,nn)
        call zfft_init_3d(nl,nm,nn,fft_struct_cz,.true.)
        nbox_cz=id*nm*nn*2
        nfft_cz=nl*nm*nn*2
        return
      entry wd_init(id,nl,nm,nn)
        call dfft_init_3d(nl,nm,nn,fft_struct_wd,.true.)
        nbox_wd=id*nm*nn
        nfft_wd=nl*nm*nn
        return
      entry wz_init(id,nl,nm,nn)
        call zfft_init_3d(nl,nm,nn,fft_struct_wz,.true.)
        nbox_wz=id*nm*nn*2
        nfft_wz=nl*nm*nn*2
        return
      end
      subroutine mulfac(afft,nbox,nfft)
      implicit real*8(a-h,o-z)
      dimension afft(nbox)
      scale=dble(nfft)
      do i=1,nbox
        afft(i)=afft(i)*scale
      enddo
      return
      end
