!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Shift_Kpoint                                       #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
   subroutine Shift_Kpoint
      use m_Kpoints,               only : vkxyz_ek, kv3_ek
      use m_Files,                 only : nfout
      use m_Epsilon_ek,            only : sw_mass, mass_dir, ikshift, mass_direction, mass_ik, norm, check_gamma_point, vkxyz_ek_org
      use m_Crystal_Structure,     only : altv
      use m_Control_Parameters,    only : printable
      use m_Const_parameters,      only : DP, CARTS, BUCS, PAI
      implicit none
      integer                         :: ik, ind
      real(DP),dimension(3)           :: ashift, rshift
      real(DP)                        :: len_ashift
      if(sw_mass ==0.or.ikshift==0.0d0) return
! backup vkxyz_ek
      vkxyz_ek_org(1:kv3_ek,1:3,1:2) = vkxyz_ek(1:kv3_ek,1:3,1:2)
! set k-point shift
      ashift(1:3) = mass_dir(1:3)
      call norm(ashift,len_ashift)
      if(printable) write(nfout,'(1x,"!* len_mass_dir =",f10.5)') len_ashift
      ashift(1:3)= ashift(1:3)*ikshift
      if(printable) write(nfout,'(1x,"!* ikshift = ",f10.5)') ikshift
      if(printable) write(nfout,'(1x,"!* ashift = (",3f10.5,")")') ashift(1:3)
      if(mass_direction==0) then
! mass_ik point shift
         ik = mass_ik
         if(printable) write(nfout,'(1x,"!* shift of ik = ",i4," for mass calculation")') ik
         call shift_ikpoint
      else
! gamma point shift
         if(printable) write(nfout,'(1x,"!* shift of gamma point for mass calculation")')
         do ik = 1, kv3_ek
            call check_gamma_point(vkxyz_ek(ik,1:3,CARTS),ind)
            if(ind/=0) exit
         end do
         call shift_ikpoint
      end if
      if(printable) write(nfout,'(1x,"!* rshift = (",3f10.5,")")') rshift(1:3)
      if(printable) write(nfout,'(1x,"!* shifted ik = ",i4," point = (",3f10.5,")")') &
     & ik, vkxyz_ek(ik,1,BUCS),vkxyz_ek(ik,2,BUCS),vkxyz_ek(ik,3,BUCS)
      contains
       subroutine shift_ikpoint
         implicit none
         real(DP), dimension(3,3) :: work
         work(1:3,1:3) = altv(1:3,1:3)/(2.0d0*PAI)
         rshift(1) = work(1,1)*ashift(1)+work(1,2)*ashift(2)+work(1,3)*ashift(3)
         rshift(2) = work(2,1)*ashift(1)+work(2,2)*ashift(2)+work(2,3)*ashift(3)
         rshift(3) = work(3,1)*ashift(1)+work(3,2)*ashift(2)+work(3,3)*ashift(3)
         vkxyz_ek(ik,1:3,BUCS) = vkxyz_ek(ik,1:3,BUCS) + rshift(1:3)
       end subroutine shift_ikpoint
   end subroutine Shift_Kpoint
