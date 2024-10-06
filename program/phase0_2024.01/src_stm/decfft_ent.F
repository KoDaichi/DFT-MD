!=======================================================================
!
!  SOFTWARE NAME : PHASE version 2.0
!
!  SUBROUINE: dxml_fft, mulfac, dummy_fft
!
!  AUTHOR(S): H. Katagiri   August/20/2003
!  
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!     The names of the contributors to the original version are as
!  follows. Dr. Hideki Katagiri, DR. Koichi Kato, Dr. Tsuyoshi Miyazaki,
!  Dr. Yoshitada Morikawa, Dr. Hideaki Sawada, Dr. Toshihiro Uchiyama,
!  Dr. Tsuyoshi Uda and Dr. Takahiro Yamasaki.
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation 
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan. 
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!
#ifdef DECFFT
      subroutine dxml_fft(iop,id,nl,nm,nn,afft)
c $Id: decfft_ent.F,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit real*8(a-h,o-z)
      include '/usr/include/DXMLDEF.FOR'
      real*8 afft(*)
      record /DXML_d_fft_structure_3D/ fft_struct_cd,fft_struct_wd
     &     , fft_struct_pwd
      record /DXML_z_fft_structure_3D/ fft_struct_cz,fft_struct_wz
     &     , fft_struct_pwz
      save fft_struct_cd,fft_struct_wd
      save fft_struct_cz,fft_struct_wz
      save nfft_cd,nfft_cz,nfft_wd,nfft_wz
      save nbox_cd,nbox_cz,nbox_wd,nbox_wz
c ------- Positron start 
      save fft_struct_pwd,fft_struct_pwz
      save nfft_pwd,nfft_pwz
      save nbox_pwd,nbox_pwz
c ------- Positron end
c
c     IOP=1..DIRECT, IOP=2..INVERSE (REAL TRANSFORM)
c     IOP=1..INVERSE, IOP=2..DIRECT (COMPLEX TRANSFORM)
c
      entry wd_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call dfft_apply_3d('r','c','f',afft,afft,id,nm,
     $                        fft_struct_wd,1,1,1)
        else
          call dfft_apply_3d('c','r','b',afft,afft,id,nm,
     $                        fft_struct_wd,1,1,1)
          call mulfac(afft,nbox_wd,nfft_wd)
        endif
      return

      entry pwd_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call dfft_apply_3d('r','c','f',afft,afft,id,nm,
     $                        fft_struct_pwd,1,1,1)
        else
          call dfft_apply_3d('c','r','b',afft,afft,id,nm,
     $                        fft_struct_pwd,1,1,1)
          call mulfac(afft,nbox_pwd,nfft_pwd)
        endif
      return
c
      entry wz_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call zfft_apply_3d('c','c','f',afft,afft,id,nm,
     $                        fft_struct_wz,1,1,1)
        else
          call zfft_apply_3d('c','c','b',afft,afft,id,nm,
     $                        fft_struct_wz,1,1,1)
          call mulfac(afft,nbox_wz,nfft_wz)
        endif
      return
c
      entry pwz_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call zfft_apply_3d('c','c','f',afft,afft,id,nm,
     $                        fft_struct_pwz,1,1,1)
        else
          call zfft_apply_3d('c','c','b',afft,afft,id,nm,
     $                        fft_struct_pwz,1,1,1)
          call mulfac(afft,nbox_pwz,nfft_pwz)
        endif
      return
c
      entry cd_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call dfft_apply_3d('r','c','f',afft,afft,id,nm,
     $                        fft_struct_cd,1,1,1)
        else
          call dfft_apply_3d('c','r','b',afft,afft,id,nm,
     $                        fft_struct_cd,1,1,1)
          call mulfac(afft,nbox_cd,nfft_cd)
        endif
      return
c
      entry cz_fft(iop,id,nm,afft)
        if(iop.eq.1) then
          call zfft_apply_3d('c','c','b',afft,afft,id,nm,
     $                        fft_struct_cz,1,1,1)
          call mulfac(afft,nbox_cz,nfft_cz)
        else
          call zfft_apply_3d('c','c','f',afft,afft,id,nm,
     $                        fft_struct_cz,1,1,1)
        endif
      return
          
      entry wd_init(id,nl,nm,nn)
        call dfft_init_3d(nl,nm,nn,fft_struct_wd,.true.)
        nbox_wd=id*nm*nn
        nfft_wd=nl*nm*nn
        return
      entry pwd_init(id,nl,nm,nn)
        call dfft_init_3d(nl,nm,nn,fft_struct_pwd,.true.)
        nbox_pwd = id*nm*nn
        nfft_pwd = nl*nm*nn
        return
      entry wz_init(id,nl,nm,nn)
        call zfft_init_3d(nl,nm,nn,fft_struct_wz,.true.)
        nbox_wz=id*nm*nn*2
        nfft_wz=nl*nm*nn
        return
      entry pwz_init(id,nl,nm,nn)
        call zfft_init_3d(nl,nm,nn,fft_structu_pwz,.true.)
        nbox_pwz = id*nm*nn*2
        nfft_pwz = nl*nm*nn
        return
      entry cd_init(id,nl,nm,nn)
        call dfft_init_3d(nl,nm,nn,fft_struct_cd,.true.)
        nbox_cd=id*nm*nn
        nfft_cd=nl*nm*nn
        return
      entry cz_init(id,nl,nm,nn)
        call zfft_init_3d(nl,nm,nn,fft_struct_cz,.true.)
        nbox_cz=id*nm*nn*2
        nfft_cz=nl*nm*nn
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
#else
      subroutine dummy_fft()
c $Id: decfft_ent.F,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      return
      end
#endif
