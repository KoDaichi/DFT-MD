!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 579 $)
!
!  MODULE: m_PAW_Hartree
!
!  AUTHOR(S): T. Yamasaki, T. Yamamoto and T. Ohno November/2009
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!#========================================================================
!
! Bug fix:  2016/10/13
!    *Antiferromagnetic calculation
!        The variable flg_done is neglected for the moment
!        because the skipping mechanism does not work properly in some cases.
!
!#========================================================================
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
module m_PAW_Hartree
  use m_Const_Parameters,     only   : DP, chg_symm_level1
  use m_Control_Parameters,   only   : nspin,af,printable,ipripp, charge_symm_mode
  use m_Ionic_System,         only   : ntyp,ivan,iatomn,iatom,natm, iwei, ityp
  use m_Charge_Density,       only   : hsr,hsro
  use m_PseudoPotential,      only   : dion_hartree,dion_hartree_now, ilmt,ipaw &
       , n_cijkclmk,ilmt3_cijkclmk,ilmt4_cijkclmk &
       , CijkClmkVVVVijlm_k &
       , CijkClmkVVVVijlm_k_ae &
       , flg_symmtry,ia2ia_symmtry_op &
       , n_cijkclmn,ilmt3_cijkclmn,ilmt4_cijkclmn &
       , CijkClmnVVVVijlm_kn &
       , CijkClmnVVVVijlm_kn_ae
  use m_Crystal_Structure,    only   : nopr
  use m_Files,                only   : nfout
  
! =========================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,   only : noncol
! ============================================================ 11.0

  implicit none
  private
  
  real(DP),pointer,dimension(:):: eh_paw,eho_paw
  
  public:: m_PAW_Hartree_alloc_ehs
  public:: m_PAWH_get_dion_hartree
  public:: m_PAWH_get_dion_hartree_now
  public:: m_PAWH_get_energy_hartree
    
contains
    
  subroutine m_PAW_Hartree_alloc_ehs
    allocate(eh_paw(natm));eh_paw=0.d0
    allocate(eho_paw(natm));eho_paw=0.d0
    return
  end subroutine m_PAW_Hartree_alloc_ehs
  
  subroutine m_PAW_Hartree_dealloc_ehs
    deallocate(eh_paw,eho_paw)
    return
  end subroutine m_PAW_Hartree_dealloc_ehs
  
  subroutine m_PAWH_get_dion_hartree(nfout)
    integer,intent(in):: nfout
    integer:: ia,ja,it,lmt1,lmt2,lmt3,lmt4,is,mp,nlmt
    real(DP):: fac,sum,sum_ae
    real(DP),pointer:: dion_hartree_ae(:,:,:)
    integer:: iopr, nopr_max
    integer, parameter :: PRINTLEVEL = 2
    logical,pointer:: flg_done(:)
    
! ========================== added by K. Tagami ================== 11.0
    integer :: ismax
    if ( noncol ) then
       ismax = 1
    else
       ismax = nspin
    endif
! ================================================================ 11.0
    
    if ( charge_symm_mode >= chg_symm_level1 ) then
       nopr_max = 1
    else
       nopr_max = nopr
    endif

    nlmt=maxval(ilmt)
    allocate(dion_hartree_ae(nlmt,nlmt,natm))
    dion_hartree_ae=0.d0
    
    do ia=1,natm
       it=ityp(ia)
!TypeLoop:   do it=1,ntyp
!                if(ityp(ia)/=it) cycle
       if(ipaw(it)/=1) then
          dion_hartree(:,:,ia)=0.d0
!                    exit TypeLoop
          cycle
       end if
       
       if(.not.flg_symmtry) then
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                sum=0.d0
                sum_ae=0.d0
                do mp=1,n_cijkclmk(lmt1,lmt2,it)
                   lmt3=ilmt3_cijkclmk(lmt1,lmt2,mp,it)
                   lmt4=ilmt4_cijkclmk(lmt1,lmt2,mp,it)
                   if(lmt3.gt.lmt4) then
                      write(nfout,*) 'Error in calcDionHartree !'
                      write(nfout,'(a,i4)') 'lmt3 = ',lmt3
                      write(nfout,'(a,i4)') 'lmt4 = ',lmt4
                      stop
                   end if
                   fac=2.d0;if(lmt3.eq.lmt4) fac=1.d0

! ========================== modified by K. Tagami ================== 11.0
!                   do is=1,nspin
                   do is=1, ismax
! =================================================================== 11.0
                      sum=sum+CijkClmkVVVVijlm_k(lmt1,lmt2,mp,it)* &
                           hsr(ia,lmt3,lmt4,is)*fac
                      sum_ae=sum_ae+CijkClmkVVVVijlm_k_ae(lmt1,lmt2,mp,it)* &
                           hsr(ia,lmt3,lmt4,is)*fac
                   end do

                end do
                dion_hartree(lmt1,lmt2,ia)=sum
                dion_hartree_ae(lmt1,lmt2,ia)=sum_ae
             end do
          end do
       else
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                sum=0.d0
                sum_ae=0.d0
                do iopr=1,nopr_max
                   do mp=1,n_cijkclmn(lmt1,lmt2,ia,iopr)
                      lmt3=ilmt3_cijkclmn(lmt1,lmt2,mp,ia,iopr)
                      lmt4=ilmt4_cijkclmn(lmt1,lmt2,mp,ia,iopr)
                      if(lmt3.gt.lmt4) then
                         write(nfout,*) 'Error in calcDionHartree !'
                         write(nfout,'(a,i4)') 'lmt3 = ',lmt3
                         write(nfout,'(a,i4)') 'lmt4 = ',lmt4
                         stop
                      end if
                      fac=2.d0;if(lmt3.eq.lmt4) fac=1.d0

! ========================== modified by K. Tagami ================== 11.0
!                      do is=1,nspin
                      do is=1, ismax
! =================================================================== 11.0

                         sum=sum+CijkClmnVVVVijlm_kn(lmt1,lmt2,mp,ia,iopr)* &
                              hsr(abs(ia2ia_symmtry_op(ia,iopr)),lmt3,lmt4,is)*fac
                         sum_ae=sum_ae+CijkClmnVVVVijlm_kn_ae(lmt1,lmt2,mp,ia,iopr)* &
                              hsr(abs(ia2ia_symmtry_op(ia,iopr)),lmt3,lmt4,is)*fac
                         !print *,CijkClmnVVVVijlm_kn_ae(lmt1,lmt2,mp,ia,iopr)
                      end do
                   end do
                end do
                dion_hartree(lmt1,lmt2,ia)=sum/dble(nopr_max)
                dion_hartree_ae(lmt1,lmt2,ia)=sum_ae/dble(nopr_max)
             end do
          end do
       end if
            
!            end do TypeLoop
    end do

!ASMS    if(af /= 0) then
!ASMS       allocate(flg_done(natm))
!ASMS       flg_done=.false.
!ASMS       do ia=1,natm
!ASMS          if(flg_done(ia)) cycle
!ASMS          ja=ia2ia_symmtry_op(ia,nopr+af)
!ASMS          dion_hartree(:,:,ia)=(dion_hartree(:,:,ia)+dion_hartree(:,:,ja))/2.d0
!ASMS          dion_hartree(:,:,ja)=dion_hartree(:,:,ia)
!ASMS          flg_done(ja)=.true.
!ASMS       end do
!ASMS       deallocate(flg_done)
!ASMS    end if
    
    if(ipripp>=PRINTLEVEL.and.printable)then
       write(nfout,*)
       write(nfout,*) ' -- dion_hartree ---'
       do ia=1,natm
          write(nfout,'(a,i2)') 'ia = ',ia
          it=ityp(ia)
          if(ipaw(it)/=1) cycle
          do lmt1 = 1, ilmt(it)
             write(nfout,'(i3,15f15.9/15f15.9)') lmt1 &
                  &               ,(dion_hartree(lmt1,lmt2,ia),lmt2 = 1, ilmt(it))
          enddo
       end do
       write(nfout,*)
       write(nfout,*) ' -- dion_hartree_ae ---'
       do ia=1,natm
          write(nfout,'(a,i2)') 'ia = ',ia
          it=ityp(ia)
          if(ipaw(it)/=1) cycle
          do lmt1 = 1, ilmt(it)
             write(nfout,'(i3,15f15.9/15f15.9)') lmt1 &
                  &               ,(dion_hartree_ae(lmt1,lmt2,ia),lmt2 = 1, ilmt(it))
          enddo
       end do
    endif
    
    deallocate(dion_hartree_ae)
!    do lmt1=1,8
!    print '(8f19.6)',(dion_hartree(lmt1,lmt2,1),lmt2=1,8)
!    end do
!    stop
    return
  end subroutine m_PAWH_get_dion_hartree
    
  subroutine m_PAWH_get_dion_hartree_now(nfout)
    integer,intent(in):: nfout
    integer:: ia,ja,it,lmt1,lmt2,lmt3,lmt4,is,mp
    real(DP):: fac,sum
    integer:: iopr, nopr_max
    logical,pointer:: flg_done(:)
    
! ========================== added by K. Tagami ================== 11.0
    integer :: ismax
    if ( noncol ) then
       ismax = 1
    else
       ismax = nspin
    endif
! ================================================================ 11.0

    if ( charge_symm_mode >= chg_symm_level1 ) then
       nopr_max = 1
    else
       nopr_max = nopr
    endif

    do ia=1,natm
       it=ityp(ia)
       !TypeLoop:   do it=1,ntyp
       !                if(ityp(ia)/=it) cycle
       if(ipaw(it)/=1) then
          dion_hartree_now(:,:,ia)=0.d0
!                    exit TypeLoop
          cycle
       end if
       
       if(.not.flg_symmtry) then
          
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                sum=0.d0
                do mp=1,n_cijkclmk(lmt1,lmt2,it)
                   lmt3=ilmt3_cijkclmk(lmt1,lmt2,mp,it)
                   lmt4=ilmt4_cijkclmk(lmt1,lmt2,mp,it)
                   if(lmt3.gt.lmt4) then
                      write(nfout,*) 'Error in calcDionHartree !'
                      write(nfout,'(a,i4)') 'lmt3 = ',lmt3
                      write(nfout,'(a,i4)') 'lmt4 = ',lmt4
                      stop
                   end if
                   fac=2.d0;if(lmt3.eq.lmt4) fac=1.d0
! ========================== modified by K. Tagami ================== 11.0
!                   do is=1,nspin
                   do is=1, ismax
! =================================================================== 11.0
                      sum=sum+CijkClmkVVVVijlm_k(lmt1,lmt2,mp,it)* &
                           hsr(ia,lmt3,lmt4,is)*fac
                   end do
                end do
                dion_hartree_now(lmt1,lmt2,ia)=sum
             end do
          end do
          
       else
          
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                sum=0.d0
                do iopr=1,nopr_max
                   do mp=1,n_cijkclmn(lmt1,lmt2,ia,iopr)
                      lmt3=ilmt3_cijkclmn(lmt1,lmt2,mp,ia,iopr)
                      lmt4=ilmt4_cijkclmn(lmt1,lmt2,mp,ia,iopr)
                      if(lmt3.gt.lmt4) then
                         write(nfout,*) 'Error in calcDionHartree !'
                         write(nfout,'(a,i4)') 'lmt3 = ',lmt3
                         write(nfout,'(a,i4)') 'lmt4 = ',lmt4
                         stop
                      end if
                      fac=2.d0;if(lmt3.eq.lmt4) fac=1.d0
! ========================== modified by K. Tagami ================== 11.0
!                   do is=1,nspin
                      do is=1, ismax
! =================================================================== 11.0
                         sum=sum+CijkClmnVVVVijlm_kn(lmt1,lmt2,mp,ia,iopr)* &
                              hsr(abs(ia2ia_symmtry_op(ia,iopr)),lmt3,lmt4,is)*fac
                      end do
                   end do
                end do
                dion_hartree_now(lmt1,lmt2,ia)=sum/dble(nopr_max)
             end do
          end do
          
       end if
                
!            end do TypeLoop
    end do
!        write(nfout,*)
!        write(nfout,*) ' -- dion_hartree ---'
!        do ia=1,natm
!            write(nfout,'(a,i2)') 'ia = ',ia
!            it=ityp(ia)
!            if(.not.ipaw(it)) cycle
!            do lmt1 = 1, ilmt(it)
!                write(nfout,'(i3,9f8.5/9f8.5)') lmt1 &
!                  &               ,(dion_hartree(lmt1,lmt2,ia),lmt2 = 1, ilmt(it))
!            enddo
!        end do
!    do lmt1=1,8
!    print '(8f19.6)',(dion_hartree(lmt1,lmt2,1),lmt2=1,8)
!    end do
!    stop

!ASMS    if(af /= 0) then
!ASMS       allocate(flg_done(natm))
!ASMS       flg_done=.false.
!ASMS       do ia=1,natm
!ASMS          if(flg_done(ia)) cycle
!ASMS          ja=ia2ia_symmtry_op(ia,nopr+af)
!ASMS          dion_hartree_now(:,:,ia)=(dion_hartree_now(:,:,ia)+dion_hartree_now(:,:,ja))/2.d0
!ASMS          dion_hartree_now(:,:,ja)=dion_hartree_now(:,:,ia)
!ASMS          flg_done(ja)=.true.
!ASMS       end do
!ASMS       deallocate(flg_done)
!ASMS    end if

    return
  end subroutine m_PAWH_get_dion_hartree_now
    
  subroutine m_PAWH_get_energy_hartree
    integer:: lmt1,lmt2,is,ia,it
    real(DP):: fac

! ========================== added by K. Tagami ================== 11.0
    integer :: ismax
    if ( noncol ) then
       ismax = 1
    else
       ismax = nspin
    endif
! ================================================================ 11.0

    do ia=1,natm
       eho_paw(ia)=0.d0
       do it=1,ntyp
          if(ityp(ia)/=it .or. ipaw(it)/=1) cycle
          
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                fac=2.d0*iwei(ia);if(lmt1.eq.lmt2) fac=iwei(ia)
! ========================== modified by K. Tagami ================== 11.0
!                do is=1,nspin
                do is=1, ismax
! =================================================================== 11.0
                   eho_paw(ia)=eho_paw(ia)+ &
                        fac*hsr(ia,lmt1,lmt2,is)* &
                        dion_hartree(lmt1,lmt2,ia)
                end do
             end do
          end do
          eho_paw(ia)=eho_paw(ia)/2.d0
       end do
    end do
    
    call m_PAWH_get_dion_hartree(nfout)
    
    do ia=1,natm
       eh_paw(ia)=0.d0
       do it=1,ntyp
          if(ityp(ia)/=it .or. ipaw(it)/=1) cycle
          
          do lmt1=1,ilmt(it)
             do lmt2=lmt1,ilmt(it)
                fac=2.d0*iwei(ia);if(lmt1.eq.lmt2) fac=iwei(ia)
! ========================== modified by K. Tagami ================== 11.0
!                do is=1,nspin
                do is=1, ismax
! =================================================================== 11.0
                   eh_paw(ia)=eh_paw(ia)+ &
                        fac*hsr(ia,lmt1,lmt2,is)* &
                        dion_hartree(lmt1,lmt2,ia)
                end do
             end do
          end do
          eh_paw(ia)=eh_paw(ia)/2.d0
       end do
    end do
    return
    
  end subroutine m_PAWH_get_energy_hartree
  

end module m_PAW_Hartree
