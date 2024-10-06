module m_OP_decomp_band
  use mpi
  implicit none
!  include 'mpif.h'
  integer istatus(mpi_status_size)

contains

  subroutine m_OP_wd_PROCAR
    use m_Control_Parameters, only : neg, nspin, ndim_spinor, ekmode, ndim_magmom
    use m_Control_Parameters, only : procar_save_memory_mode, sw_procar_full_bz, &
         &                           split_procar_file, num_procar_files_once, &
         &                           procar_sort_kpt, procar_sph_zaxis_to_magdir, &
         &                           procar_norm_zaxis_to_magdir, noncol
    use m_Const_Parameters,  only : zi, PAI, DP, CMPLDP, ON, OFF, BOHR, BUCS, HARTREE, &
         &                          GAMMA
    use m_Ionic_System,    only : ntyp, ityp, speciesname, pos, natm, &
         &                        lattice_system_from_m_CS_SG, iproj_group
    use m_PseudoPotential,   only : nlmt, ilmt, lmta, q, &
         &                          nlmta_phi, nlmtt_phi, qorb, m_PP_tell_iorb_lmtt, &
         &                          m_PP_tell_iorb_ia_l_m_tau, ilmt_phi, &
         &                          mtp_phi, lmta_phi, ltp_phi, taup_phi
    use m_PseudoPotential, only : iproj_phi, lmtt_phi, nrorb, irorb, crorb, nloc, &
         &                        nlmta_phi
    use m_Crystal_Structure, only : Global_Quantz_Axis_Fixed, sw_fix_global_quantz_axis
    use m_Crystal_structure,  only : op, nopr, altv
    use m_CD_Mag_Moment,    only : RhoMag_on_atom
    use m_Kpoints,  only : kv3, kv3_fbz, star_of_k, num_star_of_k, iopr_k_fbz_to_ibz, &
         &                 trev_k_fbz_to_ibz, vkxyz_fbz, m_Kp_set_star_of_k, &
         &                 vkxyz, k_symmetry, qwgt, kv3_ek
    use m_CS_SpaceGroup,  only : m_CS_gen_opr_rspace_full
    use m_CS_Magnetic,   only : m_CS_gen_opr_spinor_full, magmom_dir_inversion_opr_flag,&
         &                      determinant_op
    use m_Files,  only : nfout
    use m_SpinOrbit_Potential,  only : m_SO_set_MatU_ylm_RC
    use m_SpinOrbit_Potential,  only : MatU_ylm_RC_L0,  MatU_ylm_RC_L1,  &
         &                             MatU_ylm_RC_L2,  MatU_ylm_RC_L3
    use m_Parallelization,  only : ista_snl, iend_snl, ista_k, iend_k, npes, mype, &
         &                         myrank_k, myrank_e, nrank_k, MPI_CommGroup, &
         &                         ista_e, iend_e, istep_e, map_k, map_z, map_e, &
         &                         mpi_e_world, mpi_k_world
    use m_IterationNumbers,     only : nk_in_the_process

    use m_Electronic_Structure,  only : compr_l, compi_l, eko_l, occup_l, efermi
    use m_Nonlocal_Potential,   only : norm_phig

    integer :: kfac, num_orbitals
    integer :: ind_max, lun, ierr

    real(kind=DP), allocatable :: compr(:,:,:,:), compi(:,:,:,:), norm_phig_mpi(:,:)
    real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:)
    integer, allocatable :: index_p2v(:)
    character*3, allocatable :: orb_name(:)
    complex(kind=CMPLDP), allocatable :: rotated_spinor(:,:,:)
    real(kind=DP), allocatable :: rotated_sph(:,:,:,:)

    integer :: nr, num, is, i, j
    character*4 char1
    character*72 file0, file1

    integer :: lmax, mmax

    call write_poscar

    kfac = 1;    num_orbitals = nlmta_phi *kfac
!!
    if ( noncol ) then
       if ( procar_sph_zaxis_to_magdir == ON ) then
          lmax = nloc;   mmax = 2*lmax -1
          allocate( rotated_sph(natm,lmax,mmax,mmax) ); rotated_sph = 0.0d0
          call set_rotated_sph
       else if ( procar_norm_zaxis_to_magdir == ON ) then
          lmax = 2;   mmax = 2*lmax -1
          allocate( rotated_sph(natm,lmax,mmax,mmax) ); rotated_sph = 0.0d0
          call set_rotated_sph
       endif
    endif

    if ( sw_procar_full_bz == ON ) call m_Kp_set_star_of_k

    if ( split_procar_file == ON ) then
       allocate( norm_phig_mpi(nlmtt_phi,ista_snl:iend_snl) );  norm_phig_mpi = 0.0d0
       allocate( compr( neg, num_orbitals,1, ista_k:iend_k ) ); compr = 0.0d0
       allocate( compi( neg, num_orbitals,1, ista_k:iend_k ) ); compi = 0.0d0
       call set_compri_etc_split
       if (noncol) then
          allocate( eko_wk(neg,ista_snl:iend_snl) );  eko_wk = 0.0d0
          allocate( occ_wk(neg,ista_snl:iend_snl) );  occ_wk = 0.0d0
       else
          allocate( eko_wk(neg,ista_k:iend_k) );  eko_wk = 0.0d0
          allocate( occ_wk(neg,ista_k:iend_k) );  occ_wk = 0.0d0
       endif
       call set_eko_occup_wk_split

    else
       allocate( norm_phig_mpi(nlmtt_phi,kv3/nspin) );  norm_phig_mpi = 0.0d0
       select case (procar_save_memory_mode)
       case (0)
          allocate( compr( neg, num_orbitals,1, kv3 ) ); compr = 0.0d0
          allocate( compi( neg, num_orbitals,1, kv3 ) ); compi = 0.0d0
          call set_compri_etc
       case (1)
          if ( mype == 0 ) then
             allocate( compr( neg, num_orbitals,1, kv3 ) ); compr = 0.0d0
             allocate( compi( neg, num_orbitals,1, kv3 ) ); compi = 0.0d0
          endif
          call set_compri_etc2
       end select

       allocate( eko_wk(neg,kv3/ndim_spinor) );  eko_wk = 0.0d0
       allocate( occ_wk(neg,kv3/ndim_spinor) );  occ_wk = 0.0d0
       call set_eko_occup_wk
    endif

    ind_max = 9          ! up to d-orbital
    allocate( orb_name(ind_max) );   allocate( index_p2v(ind_max) )
    call set_orb_name(ind_max,orb_name)
    call set_index_converter_p2v(ind_max,index_p2v)

    if ( split_procar_file == ON ) then
       nr = nrank_k /num_procar_files_once +1
       num = 0

       lun = 10000 +myrank_k
       write(char1,'(I4.4)') myrank_k

       iloop: Do i=1, nr
          Do j=1, num_procar_files_once
             num = num +1
             if ( num > npes )    goto 200
             if ( num /= myrank_k +1 ) cycle
             if ( myrank_e /= 0 )   cycle

             if ( noncol ) then
                allocate( rotated_spinor(natm,ndim_spinor,ndim_spinor) )
                call set_spinor_basis

                file0 = "PROCAR_phase"
                file1 = trim(adjustl(file0)) // '.' // char1
                open( lun, file=file1, status="unknown", form="formatted" )
                if ( sw_procar_full_bz == ON ) then
                   call write_file_noncl_fbz(lun)
                else
                   call write_file_noncl(lun)
                endif
                close(lun);   deallocate( rotated_spinor )

             else
                Do is=1, nspin
                   if ( nspin == 1 ) then
                      file0 = "PROCAR_phase"
                   else
                      if ( is == 1 ) then
                         file0 = "PROCAR_phase.up"
                      else
                         file0 = "PROCAR_phase.dn"
                      endif
                   endif

                   file1 = trim(adjustl(file0)) // '.' // char1
                   open( lun, file=file1, status="unknown", form="formatted" )
                   if ( sw_procar_full_bz == ON ) then
                      call write_file_col_fbz(lun,is);
                   else
                      call write_file_col(lun,is);
                   endif
                   close(lun)
                End Do
             endif
          End Do
200       continue
          call mpi_barrier( MPI_CommGroup, ierr )
       End Do iloop

    else if ( mype == 0 ) then
       lun = 1000
#if 0
       open( lun, file="PROCAR_phase", status="unknown", form="formatted")
#else
       if ( ekmode == OFF ) then
          open( lun, file="PROCAR_phase", status="unknown", form="formatted")
       else
          if ( nk_in_the_process == 1 ) then
             open( lun, file="PROCAR_phase", status="unknown", form="formatted")
          else
             open( lun, file="PROCAR_phase", status="unknown", form="formatted", &
                  &     position="append" )
          endif
       endif
#endif
       if ( noncol ) then
          allocate( rotated_spinor(natm,ndim_spinor,ndim_spinor) )
          call set_spinor_basis

          if ( sw_procar_full_bz == ON ) then
             if ( procar_sort_kpt == ON ) then
                call write_file_noncl_fbz_sort(lun)
             else
                call write_file_noncl_fbz(lun)
             endif
          else
             call write_file_noncl(lun)
          endif
          deallocate( rotated_spinor )

       else
          Do is=1, nspin
             if ( sw_procar_full_bz == ON ) then
                if ( procar_sort_kpt == ON ) then
                   call write_file_col_fbz_sort(lun,is)
                else
                   call write_file_col_fbz(lun,is)
                endif
             else
                call write_file_col(lun,is)
             endif
          End Do
       endif
       close(lun)

    endif

    deallocate(orb_name);  deallocate(index_p2v)
    if ( allocated(compr) ) deallocate( compr );
    if ( allocated(compi) ) deallocate( compi )
    deallocate( norm_phig_mpi )

    deallocate( eko_wk );    deallocate( occ_wk )
    if ( allocated(rotated_sph) ) deallocate( rotated_sph )

  contains

    subroutine calc_euler_angles2( axis, alpha, beta, gamma )   ! Z-Y-Z
      real(kind=DP), intent(in) :: axis(3,3)
      real(kind=DP), intent(out) :: alpha, beta, gamma
      !
      real(kind=DP) :: sigma, delta, c1, c2, c3, c4, c5
      real(kind=DP) :: mat(3,3)
      !
      mat = axis
      if ( mat(3,3) < 1.0D0 ) then
         if ( mat(3,3) > -1.0D0 ) then
            beta = acos( mat(3,3) )
            alpha = atan2( mat(2,3), mat(1,3) )
            gamma = atan2( mat(3,2), -mat(3,1) )
         else
            beta = PAI
            alpha = atan2( mat(2,1), mat(2,2) )
            gamma = 0.0D0;
         endif
      else
         beta = 0.0D0
         alpha = atan2( mat(2,1), mat(2,2) )
         gamma = 0.0D0
      endif
      !
      mat(1,1) = cos(alpha)*cos(beta)*cos(gamma) -sin(alpha)*sin(gamma)
      mat(1,2) = -cos(alpha)*cos(beta)*sin(gamma) -sin(alpha)*cos(gamma)
      mat(1,3) = cos(alpha)*sin(beta)
      mat(2,1) = sin(alpha)*cos(beta)*cos(gamma) +cos(alpha)*sin(gamma)
      mat(2,2) = -sin(alpha)*cos(beta)*sin(gamma) +cos(alpha)*cos(gamma)
      mat(2,3) = sin(alpha)*sin(beta)
      mat(3,1) = -sin(beta)*cos(gamma)
      mat(3,2) = sin(beta)*sin(gamma)
      mat(3,3) = cos(beta)
    end subroutine calc_euler_angles2

    subroutine set_rotated_sph  ! noncol
      integer :: ia, il1, im1, im2, im3, im4, val_l
      real(kind=DP) :: mx, my, mz, mnorm, c1
      real(kind=DP) :: theta, phi
      real(kind=DP) :: spn_quant_dir(3), angle(3), axis(3,3)
      complex(kind=CMPLDP) :: z1, z2, zsum
      complex(kind=CMPLDP), allocatable :: mat(:,:)

      integer :: lun
      character*72 file1

      call m_SO_set_MatU_ylm_RC
      if ( mype == 0 ) then
         lun = 1000;  file1 = './PROCAR_Rotated_Axis'
         open( lun, file=file1, status="unknown", form="formatted" )
         write(lun,'(A,4X,A,14X,A,14X,A)') " AtomID",  "X-axis", "Y-axis", "Z-axis"
      endif

      Do ia=1, natm
         if ( sw_fix_global_quantz_axis == ON ) then
            spn_quant_dir(1:3) = Global_Quantz_Axis_Fixed(1:3)
         else
            mx = RhoMag_on_atom(ia,2)
            my = RhoMag_on_atom(ia,3)
            mz = RhoMag_on_atom(ia,4)
            mnorm = sqrt( mx**2 +my**2 +mz**2 )
            if ( mnorm < 1.0D-8 ) then
               spn_quant_dir(1:2) = 0.0d0; spn_quant_dir(3) = 1.0d0
            else
               spn_quant_dir(1) = mx /mnorm
               spn_quant_dir(2) = my /mnorm
               spn_quant_dir(3) = mz /mnorm
            endif
         end if
!
         call set_rotated_axis( spn_quant_dir, axis )
         call calc_euler_angles2( axis, angle(1), angle(2), angle(3) )

         if ( mype == 0 ) then
            write(lun,'(I5)') ia
            write(lun,'(3F20.10)') axis(1,1:3)
            write(lun,'(3F20.10)') axis(2,1:3)
            write(lun,'(3F20.10)') axis(3,1:3)
         endif

         Do il1=1, lmax
            val_l = il1 -1
            allocate( mat(-val_l:val_l,-val_l:val_l) );  mat = 0.0d0

            Do im1=-val_l ,val_l
               Do im2=-val_l, val_l
                  call calc_small_Wigner_function( dble(val_l), dble(im1), dble(im2), &
                       &                           -angle(2), c1 )
                  z1 = exp( zi *im1 *angle(1) )
                  z2 = exp( zi *im2 *angle(3) )
                  mat(im1,im2) = z1 *c1 *z2
               End do
            End do

            Do im1=1, 2*val_l+1
               Do im2=1, 2*val_l +1
                  zsum = 0.0d0

                  select case (val_l)
                  case (0)
                     Do im3=-val_l, val_l
                        Do im4=-val_l, val_l
                           zsum = zsum +MatU_ylm_RC_L0(im1,im3)*mat(im3,im4) &
                                &       *conjg(MatU_ylm_RC_L0(im2,im4))
                        End Do
                     End Do
                  case (1)
                     Do im3=-val_l, val_l
                        Do im4=-val_l, val_l
                           zsum = zsum +MatU_ylm_RC_L1(im1,im3)*mat(im3,im4) &
                                &       *conjg(MatU_ylm_RC_L1(im2,im4))
                        End Do
                     End Do
                  case (2)
                     Do im3=-val_l, val_l
                        Do im4=-val_l, val_l
                           zsum = zsum +MatU_ylm_RC_L2(im1,im3)*mat(im3,im4) &
                                &       *conjg(MatU_ylm_RC_L2(im2,im4))
                        End Do
                     End Do
                  case (3)
                     Do im3=-val_l, val_l
                        Do im4=-val_l, val_l
                           zsum = zsum +MatU_ylm_RC_L3(im1,im3)*mat(im3,im4) &
                                &       *conjg(MatU_ylm_RC_L3(im2,im4))
                        End Do
                     End Do
                  end select
                  rotated_sph(ia,il1,im1,im2) = zsum
               End Do
            End Do
            deallocate( mat )
         End Do
      End Do
      if ( mype == 0 ) close(lun)
#if 1
      if ( mype == 0 ) then
         lun = 1000;  file1 = './PROCAR_Rotated_Sph'
         open( lun, file=file1, status="unknown", form="formatted" )
         Do ia=1, natm
            write(lun,*) "Atom ID: ", ia
            Do il1=1, lmax
               val_l = il1 -1
               write(lun,*) "l = ", val_l
               Do im1=1, 2*val_l +1
                  Do im2=1, 2*val_l +1
                     write(lun,'(F10.4)',advance='no') rotated_sph(ia,il1,im1,im2)
                  End Do
                  write(lun,*)
               End Do
            End Do
         End Do
         close(lun)
      endif
#endif
    end subroutine set_rotated_sph

    subroutine set_rotated_axis( moment, axis )
      real(kind=DP), intent(in) :: moment(3)
      real(kind=DP), intent(out) :: axis(3,3)
!
      real(kind=DP) :: norm_xy, dx, dy, dz, th, cos_th, sin_th, c1
      real(kind=DP) :: rot_mat(3,3)
      real(kind=DP), parameter :: limit = 1.0D-6
!
! Rodrigues
      norm_xy = sqrt( moment(1)**2 +moment(2)**2 )
      if ( norm_xy < limit ) then
         axis = 0.0d0
         axis(1,1) = 1.0d0;   axis(2,2) = 1.0d0;    axis(3,3) = 1.0d0
         return
      endif
!
      dx = -moment(2)/norm_xy
      dy =  moment(1)/norm_xy
      dz =  0.0d0
!
      th = atan2( norm_xy, moment(3) )
      cos_th = cos(th);   sin_th = sin(th)
      !
      c1 = 1.0d0 -cos_th
!
!!      rot_mat(1,1) = dz*dz *c1 +cos_th
      rot_mat(1,1) = dx*dx *c1 +cos_th
      rot_mat(1,2) = dx*dy *c1 -dz *sin_th
      rot_mat(1,3) = dx*dz *c1 +dy *sin_th

      rot_mat(2,1) = dy*dx *c1 +dz *sin_th
      rot_mat(2,2) = dy*dy *c1 +cos_th
      rot_mat(2,3) = dy*dz *c1 -dx *sin_th

      rot_mat(3,1) = dz*dx *c1 -dy *sin_th
      rot_mat(3,2) = dz*dy *c1 +dx *sin_th
      rot_mat(3,3) = dz*dz *c1 +cos_th
!
      axis(1:3,1) = rot_mat(1:3,1)
      axis(1:3,2) = rot_mat(1:3,2)
      axis(1:3,3) = rot_mat(1:3,3)
    end subroutine set_rotated_axis

    subroutine write_poscar
      integer :: it, ia, num

      lun = 1000
      open( lun, file="POSCAR_phase", status="unknown", form="formatted")

      write(lun,'(X,A)') "PROCAR generated by phase"
      write(lun,'(3X,F20.16)') 1.0d0
      write(lun,'(5X,3F22.16)') altv(1:3,1) *Bohr
      write(lun,'(5X,3F22.16)') altv(1:3,2) *Bohr
      write(lun,'(5X,3F22.16)') altv(1:3,3) *Bohr
      Do it=1, ntyp
         write(lun,'(3X,A)',advance='no') speciesname(it)
      End Do
      write(lun,*)
      Do it=1, ntyp
         num = 0
         Do ia=1, natm
            if ( ityp(ia) == it ) num = num +1
         End Do
         write(lun,'(3X,I4)',advance='no') num
      End Do
      write(lun,*)

      write(lun,'(A)') "Direct"
      Do it=1, ntyp
         num = 0
         Do ia=1, natm
            if ( ityp(ia) /= it ) cycle
            write(lun,'(3F22.16)') pos(ia,1:3)
         End Do
      End Do
      write(lun,*)
      Do ia=1, natm
         write(lun,'(3F22.16)') 0.0d0, 0.0d0, 0.0d0
      End Do
      close(lun)
    end subroutine write_poscar

    subroutine set_compri_etc
      integer :: ik, iksnl, ie, ib, ierr
      real(kind=DP), allocatable :: compr_mpi(:,:,:,:), compi_mpi(:,:,:,:)
      real(kind=DP), allocatable :: norm_phig_mpi2(:,:)

      integer :: ia, it, il1, il2, im1, im2, tau1, tau2, iorb1, iorb2, lmt1, lmt2
      real(kind=DP) :: ctmp_r, ctmp_i
      real(kind=DP), allocatable :: cwork_r(:,:,:), cwork_i(:,:,:)

      if ( noncol .and. procar_sph_zaxis_to_magdir == ON ) then
         allocate( cwork_r(nlmta_phi,1,1) );
         allocate( cwork_i(nlmta_phi,1,1) );

         do ik = 1, kv3
            if(map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin + 1
            do ie = ista_e, iend_e, istep_e
               ib = map_z(ie)
               do iorb1=1, nlmta_phi
                  cwork_r(iorb1,1,1) = compr_l(ib,iorb1,1,ik)
                  cwork_i(iorb1,1,1) = compi_l(ib,iorb1,1,ik)
               end do

               Do ia=1, natm
                  it = ityp(ia)
                  if (iproj_group(ia) == 0) cycle

                  do lmt1=1,ilmt_phi(it)
                     il1 = ltp_phi(lmt1,it);      im1 = mtp_phi(lmt1,it)
                     tau1 = taup_phi(lmt1,it);    iorb1 = lmta_phi(lmt1,ia)
!                     immax = 2*il1 -1

                     ctmp_r = 0.0d0;   ctmp_i = 0.0d0
                     do lmt2=1, ilmt_phi(it)
                        il2 = ltp_phi(lmt2,it);    im2 = mtp_phi(lmt2,it)
                        tau2 = taup_phi(lmt2,it);  iorb2 = lmta_phi(lmt2,ia)
                        if ( il1 /= il2 ) cycle
                        if ( tau1 /= tau2 ) cycle
                        ctmp_r = ctmp_r &
                             &    +rotated_sph(ia,il1,im2,im1) &
                             &                 *cwork_r(iorb2,1,1)
                        ctmp_i = ctmp_i &
                             &    +rotated_sph(ia,il1,im2,im1) &
                             &                 *cwork_i(iorb2,1,1)
                     end do
                     compr(ie,iorb1,1,ik) = ctmp_r;   compi(ie,iorb1,1,ik) = ctmp_i
                  end do
               end Do
            end do
            norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
         end do
         deallocate( cwork_r );  deallocate( cwork_i )
      else
         do ik = 1, kv3
            if(map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin + 1
            do ie = ista_e, iend_e, istep_e
               ib = map_z(ie)
               compr(ie,1:num_orbitals,1,ik) = compr_l(ib,1:num_orbitals,1,ik)
               compi(ie,1:num_orbitals,1,ik) = compi_l(ib,1:num_orbitals,1,ik)
            end do
            norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
         end do
      endif

      if ( npes >1 ) then
         allocate( compr_mpi( neg, num_orbitals, 1, kv3 ) ); compr_mpi = 0.0d0
         allocate( compi_mpi( neg, num_orbitals, 1, kv3 ) ); compi_mpi = 0.0d0
         allocate( norm_phig_mpi2( nlmtt_phi, kv3/nspin ) ); norm_phig_mpi2 = 0.0d0
         call mpi_allreduce( compr, compr_mpi, neg*num_orbitals*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( compi, compi_mpi, neg*num_orbitals*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*kv3/nspin, &
              &              mpi_double_precision, mpi_sum, mpi_e_world(myrank_e), ierr )
         compr = compr_mpi;   compi = compi_mpi
         norm_phig_mpi = norm_phig_mpi2
         deallocate( compr_mpi ); deallocate( compi_mpi );
         deallocate( norm_phig_mpi2 )
      end if
    end subroutine set_compri_etc

    subroutine set_compri_etc2    ! save_memory
      integer :: ik, iksnl, ie, ib, ierr, iorb
      integer :: ia, it, il1, im1, tau1, il2, im2, tau2, iorb2, lmt2
      real(kind=DP) :: ctmp_r, ctmp_i
      real(kind=DP), allocatable :: compr_mpi(:,:,:,:), compi_mpi(:,:,:,:)
      real(kind=DP), allocatable :: norm_phig_mpi2(:,:)
      real(kind=DP), allocatable :: wkr(:,:), wkr_mpi(:,:), wki(:,:), wki_mpi(:,:)

      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         iksnl = (ik-1)/nspin + 1
         norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
      end do
      if ( npes >1 ) then
         allocate( norm_phig_mpi2( nlmtt_phi, kv3/nspin ) ); norm_phig_mpi2 = 0.0d0
         call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*kv3/nspin, &
              &              mpi_double_precision, mpi_sum, mpi_e_world(myrank_e), ierr )
         norm_phig_mpi = norm_phig_mpi2
         deallocate( norm_phig_mpi2 )
      end if

      allocate( wkr(neg,kv3) );    allocate( wki(neg,kv3) );

      Do iorb=1, num_orbitals
         wkr = 0.0d0;  wki = 0.0d0

         if ( noncol .and. procar_sph_zaxis_to_magdir == ON ) then
            call m_PP_tell_iorb_ia_l_m_tau( iorb, ia, il1, im1, tau1 )
            it = ityp(ia)

            do ik = 1, kv3
               if(map_k(ik) /= myrank_k) cycle
               iksnl = (ik-1)/nspin + 1
               do ie = ista_e, iend_e, istep_e
                  ib = map_z(ie)

                  ctmp_r = 0.0d0;   ctmp_i = 0.0d0
                  do lmt2=1, ilmt_phi(it)
                     il2 = ltp_phi(lmt2,it);    im2 = mtp_phi(lmt2,it)
                     tau2 = taup_phi(lmt2,it);  iorb2 = lmta_phi(lmt2,ia)
                     if ( il1 /= il2 ) cycle
                     if ( tau1 /= tau2 ) cycle
                     ctmp_r = ctmp_r &
                          &    +rotated_sph(ia,il1,im2,im1) &
                          &                 *compr_l(ib,iorb2,1,ik)
                     ctmp_i = ctmp_i &
                          &    +rotated_sph(ia,il1,im2,im1) &
                          &                 *compi_l(ib,iorb2,1,ik)
                  end do
                  wkr(ie,ik) = ctmp_r;    wki(ie,ik) = ctmp_i
               end do
            end do
         else
            do ik = 1, kv3
               if(map_k(ik) /= myrank_k) cycle
               iksnl = (ik-1)/nspin + 1
               do ie = ista_e, iend_e, istep_e
                  ib = map_z(ie)
                  wkr(ie,ik) = compr_l(ib,iorb,1,ik)
                  wki(ie,ik) = compi_l(ib,iorb,1,ik)
               end do
            end do
         endif

         if ( npes >1 ) then
            allocate( wkr_mpi( neg, kv3 ) ); wkr_mpi = 0.0d0
            allocate( wki_mpi( neg, kv3 ) ); wki_mpi = 0.0d0
            call mpi_allreduce( wkr, wkr_mpi, neg*kv3, &
                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            call mpi_allreduce( wki, wki_mpi, neg*kv3, &
                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            wkr = wkr_mpi;   wki = wki_mpi
            deallocate( wkr_mpi ); deallocate( wki_mpi );
         end if
         if ( mype == 0 ) then
            compr(:,iorb,1,:) = wkr(:,:);     compi(:,iorb,1,:) = wki(:,:)
         endif
      End Do
      deallocate( wkr );    deallocate( wki )
    end subroutine set_compri_etc2

    subroutine set_compri_etc_split
      integer :: ik, iksnl, ie, ib, ierr, iorb
      integer :: ia, it, il1, im1, tau1, il2, im2, tau2, iorb2, lmt2
      real(kind=DP) :: ctmp_r, ctmp_i
      real(kind=DP), allocatable :: compr_mpi(:,:,:,:), compi_mpi(:,:,:,:)
      real(kind=DP), allocatable :: norm_phig_mpi2(:,:)
      real(kind=DP), allocatable :: wkr(:,:), wkr_mpi(:,:), wki(:,:), wki_mpi(:,:)

      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         iksnl = (ik-1)/nspin + 1
         norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
      end do

      allocate( wkr(neg,ista_k:iend_k) );    allocate( wki(neg,ista_k:iend_k) );
      Do iorb=1, num_orbitals
         wkr = 0.0d0;  wki = 0.0d0

         if ( noncol .and. procar_sph_zaxis_to_magdir == ON ) then
            call m_PP_tell_iorb_ia_l_m_tau( iorb, ia, il1, im1, tau1 )
            it = ityp(ia)

            do ik = ista_k, iend_k
               iksnl = (ik-1)/nspin + 1
               do ie = ista_e, iend_e, istep_e
                  ib = map_z(ie)

                  ctmp_r = 0.0d0;   ctmp_i = 0.0d0
                  do lmt2=1, ilmt_phi(it)
                     il2 = ltp_phi(lmt2,it);    im2 = mtp_phi(lmt2,it)
                     tau2 = taup_phi(lmt2,it);  iorb2 = lmta_phi(lmt2,ia)
                     if ( il1 /= il2 ) cycle
                     if ( tau1 /= tau2 ) cycle
                     ctmp_r = ctmp_r &
                          &    +rotated_sph(ia,il1,im2,im1) &
                          &                 *compr_l(ib,iorb2,1,ik)
                     ctmp_i = ctmp_i &
                          &    +rotated_sph(ia,il1,im2,im1) &
                          &                 *compi_l(ib,iorb2,1,ik)
                  end do
                  wkr(ie,ik) = ctmp_r;    wki(ie,ik) = ctmp_i
               end do
            end do
         else
            do ik = ista_k, iend_k
               iksnl = (ik-1)/nspin + 1
               do ie = ista_e, iend_e, istep_e
                  ib = map_z(ie)
                  wkr(ie,ik) = compr_l(ib,iorb,1,ik)
                  wki(ie,ik) = compi_l(ib,iorb,1,ik)
               end do
            end do
         endif
         if ( npes >1 ) then
            allocate( wkr_mpi( neg, ista_k:iend_k ) ); wkr_mpi = 0.0d0
            allocate( wki_mpi( neg, ista_k:iend_k ) ); wki_mpi = 0.0d0
            call mpi_allreduce( wkr, wkr_mpi, neg*(iend_k-ista_k+1), &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_k_world(myrank_k), ierr )
            call mpi_allreduce( wki, wki_mpi, neg*(iend_k-ista_k+1), &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_k_world(myrank_k), ierr )
            wkr = wkr_mpi;   wki = wki_mpi
            deallocate( wkr_mpi ); deallocate( wki_mpi );
         end if
         compr(:,iorb,1,:) = wkr(:,:);     compi(:,iorb,1,:) = wki(:,:)
      End Do
      deallocate( wkr );    deallocate( wki )
    end subroutine set_compri_etc_split

    subroutine set_eko_occup_wk
      integer :: ik, ie, ierr, iksnl
      real(kind=DP), allocatable :: eko_mpi(:,:), occ_mpi(:,:)

      if ( noncol ) then
         do ik = 1, kv3, ndim_spinor
            if (map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin + 1
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,iksnl) = eko_l(map_z(ie),ik)
               occ_wk(ie,iksnl) = occup_l(map_z(ie),ik)
            end do
         end do
      else
         do ik = 1, kv3
            if (map_k(ik) /= myrank_k) cycle
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,ik) = eko_l(map_z(ie),ik)
               occ_wk(ie,ik) = occup_l(map_z(ie),ik)
            end do
         end do
      endif

      if ( npes > 1 ) then
         allocate( eko_mpi(neg,kv3/ndim_spinor) ); eko_mpi = 0.0d0
         allocate( occ_mpi(neg,kv3/ndim_spinor) ); occ_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*kv3/ndim_spinor, mpi_double_precision,&
              &              mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( occ_wk, occ_mpi, neg*kv3/ndim_spinor, mpi_double_precision,&
              &              mpi_sum, MPI_CommGroup, ierr )
         eko_wk = eko_mpi;   occ_wk = occ_mpi
         deallocate( eko_mpi );      deallocate( occ_mpi )
      endif
    end subroutine set_eko_occup_wk

    subroutine set_eko_occup_wk_split
      integer :: ik, ie, ierr, iksnl
      integer :: istart, iend
      real(kind=DP), allocatable :: eko_mpi(:,:), occ_mpi(:,:)

      if ( noncol ) then
         do ik = 1, kv3, ndim_spinor
            if (map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin + 1
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,iksnl) = eko_l(map_z(ie),ik)
               occ_wk(ie,iksnl) = occup_l(map_z(ie),ik)
            end do
         end do
      else
         do ik = 1, kv3
            if (map_k(ik) /= myrank_k) cycle
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,ik) = eko_l(map_z(ie),ik)
               occ_wk(ie,ik) = occup_l(map_z(ie),ik)
            end do
         end do
      endif

      if ( noncol ) then
         istart = ista_snl;  iend = iend_snl
      else
         istart = ista_k;  iend = iend_k
      endif
      if ( npes > 1 ) then
         allocate( eko_mpi(neg,istart:iend) ); eko_mpi = 0.0d0
         allocate( occ_mpi(neg,istart:iend) ); occ_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend-istart+1), &
              &              mpi_double_precision,&
              &              mpi_sum, mpi_k_world(myrank_k), ierr )
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend-istart+1), &
              &              mpi_double_precision,&
              &              mpi_sum, mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi;   occ_wk = occ_mpi
         deallocate( eko_mpi );      deallocate( occ_mpi )
      endif
    end subroutine set_eko_occup_wk_split

    subroutine write_file_col( lun, is )
      integer, intent(in) :: lun, is

      integer :: ik, iksnl, ia, it, ib
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind
      integer :: ikstart, ikend
      real(kind=DP) :: c1, weight
      complex(kind=CMPLDP) :: z1
      real(kind=DP), allocatable :: norm2(:,:), softpart(:,:,:)

      if ( mype ==0 .and. is==1 ) then
         write(lun,'(A)') "PROJ  (orbital) lm decomposed + phase factor"
      endif

      allocate( norm2(natm+1,ind_max+1) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,2) );  softpart = 0.0d0

      if ( mype ==0 ) then
         write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
         write(lun,*)
      endif

      if ( split_procar_file == ON ) then
         ikstart = ista_k +is -1;    ikend = iend_k
      else
         ikstart = is;               ikend = kv3
      endif

      Do ik=ikstart, ikend, nspin
         iksnl = (ik-1)/nspin + 1
         write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
              &     " k-point", iksnl, " :", vkxyz(ik,1:3,BUCS), &
              &     " weight = ", qwgt(ik)*nspin
         write(lun,*)

         weight = kv3 *qwgt(ik) /dble(ndim_spinor)
         if ( nspin == 1 ) weight = weight /2.0d0

         Do ib=1, neg
            softpart = 0.0d0;  norm2 = 0.0d0

            write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                 &     "band ", ib, &
                 &     " # energy ", ( eko_wk(ib,ik) -efermi )*Hartree, &
                 &     " # occ.", occ_wk(ib,ik)/weight
            write(lun,*)
            Do ia=1, natm
               it = ityp(ia)
               if(iproj_group(ia) == 0) cycle

               do lmt=1,ilmt_phi(it)
                  il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                  tau = taup_phi(lmt,it)

                  ind = (il-1)**2 +im
                  iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                  z1 = dcmplx( compr(ib,iorb,1,ik), compi(ib,iorb,1,ik) )
                  c1 = z1*conjg(z1)*( 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )

                  norm2(ia,ind) = c1
                  softpart(ia,ind,1) = compr(ib,iorb,1,ik)
                  softpart(ia,ind,2) = compi(ib,iorb,1,ik)
               end do
            End Do
!
            Do ia=1, natm
               c1 = 0.0d0
               Do ind=1, ind_max
                  c1 = c1 +norm2(ia,ind)
               End Do
               norm2(ia,ind_max+1) = c1
            End Do
            !
            Do ind=1, ind_max+1
               c1 = 0.0d0
               Do ia=1, natm
                  c1 = c1 +norm2(ia,ind)
               End Do
               norm2(natm+1,ind) = c1
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*) "   tot"

            Do ia=1, natm
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind))
               End Do
               write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1)
               write(lun,*)
            End Do

            write(lun,'(A)',advance='no') "tot "
            Do ind=1, ind_max
               write(lun,'(F18.12)',advance='no') norm2(natm+1,index_p2v(ind))
            End Do
            write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1)
            write(lun,*)

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*)

            Do ia=1, natm
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),1)
               End Do
               write(lun,*)
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),2)
               End Do
               write(lun,*)
            End Do
            write(lun,*)

         End Do
         if ( iksnl /= (kv3/nspin) ) write(lun,*)
      End Do
      deallocate(norm2); deallocate(softpart)

    end subroutine write_file_col

    subroutine write_file_col_fbz( lun, is )
      integer, intent(in) :: lun, is

      integer :: ik, iksnl, ia, it, ib, jksnl, ii
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind, lmtt2
      real(kind=DP) :: c1, weight, sign1
      complex(kind=CMPLDP) :: z1
      real(kind=DP), allocatable :: norm2(:,:), softpart(:,:,:)
      real(DP),allocatable :: crotylm(:,:,:)
      integer :: ikstart, ikend, jk, iopr, trev, mm, jorb
      real(kind=DP), allocatable :: compr_wk(:), compi_wk(:)

      if ( mype ==0 .and. is==1 ) then
         write(lun,'(A)') "PROJ  (orbital) lm decomposed + phase factor"
      endif

      allocate( norm2(natm+1,ind_max+1) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,2) );  softpart = 0.0d0
      if ( mype ==0 ) then
         write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3_fbz/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
         write(lun,*)
      endif

      if ( split_procar_file == ON ) then
         ikstart = ista_k +is -1;    ikend = iend_k
      else
         ikstart = is;               ikend = kv3
      endif

      allocate( compr_wk(nlmta_phi) )
      allocate( compi_wk(nlmta_phi) )

      Do ik=ikstart, ikend, nspin
         iksnl = (ik-1)/nspin + 1

         weight = kv3 *qwgt(ik) /dble(ndim_spinor)
         if ( nspin == 1 ) weight = weight /2.0d0

         Do ii=1, num_star_of_k(ik)
            jk = star_of_k(ik,ii)
            jksnl = (jk-1)/nspin + 1
            iopr = iopr_k_fbz_to_ibz(jk)
            trev = trev_k_fbz_to_ibz(jk)

            write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
                 &     " k-point", jksnl, " :", vkxyz_fbz(jk,1:3,BUCS), &
                 &     " weight = ", 1.0d0 / (kv3_fbz/nspin)
            write(lun,*)

            Do ib=1, neg
               softpart = 0.0d0;  norm2 = 0.0d0

               write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                    &     "band ", ib, &
                    &     " # energy ", ( eko_wk(ib,ik) -efermi )*Hartree, &
                    &     " # occ.", occ_wk(ib,ik)/weight
               write(lun,*)
               !
               if ( k_symmetry(ik) == GAMMA ) then
               else
                  compr_wk = 0.d0;          compi_wk = 0.d0
                  do iorb=1,nlmta_phi
                     do mm=1,nrorb(iorb,iopr)
                        jorb = irorb(mm,iorb,iopr)
                        compr_wk(iorb) = compr_wk(iorb) +compr(ib,jorb,1,ik) &
                             &          *crorb(mm,iorb,iopr)
                        compi_wk(iorb) = compi_wk(iorb) +compi(ib,jorb,1,ik) &
                             &          *crorb(mm,iorb,iopr)
                     end do
                  end do
               endif
!
               Do ia=1, natm
                  it = ityp(ia)
                  if (iproj_group(ia) == 0) cycle

                  do lmt=1,ilmt_phi(it)
                     il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                     tau = taup_phi(lmt,it)
                     ind = (il-1)**2 +im
                     iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                     z1 = dcmplx( compr_wk(iorb), compi_wk(iorb) )

                     sign1 = 1.0d0
                     if ( trev == 1 ) sign1 = -sign1
#if 0
                     if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) sign1 = -sign1
#endif
                     if ( sign1 < .0d0 ) z1 = conjg(z1)

#if 1
                     if ( determinant_op(iopr) < 0 ) then
                        z1 = z1 *( (-1)**(il -1) )
                     endif
#endif
                     c1 = z1*conjg(z1)*( 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )

                     norm2(ia,ind) = c1
                     softpart(ia,ind,1) = real( z1 )
                     softpart(ia,ind,2) = aimag( z1 )
                  end do
               End Do
!
               Do ia=1, natm
                  c1 = 0.0d0
                  Do ind=1, ind_max
                     c1 = c1 +norm2(ia,ind)
                  End Do
                  norm2(ia,ind_max+1) = c1
               End Do
!
               Do ind=1, ind_max+1
                  c1 = 0.0d0
                  Do ia=1, natm
                     c1 = c1 +norm2(ia,ind)
                  End Do
                  norm2(natm+1,ind) = c1
               End Do

               write(lun,'(A)',advance='no') "ion "
               Do ind=1, ind_max
                  write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
               End Do
               write(lun,*) "   tot"

               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind))
                  End Do
                  write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1)
                  write(lun,*)
               End Do

               write(lun,'(A)',advance='no') "tot "
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') norm2(natm+1,index_p2v(ind))
               End Do
               write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1)
               write(lun,*)

               write(lun,'(A)',advance='no') "ion "
               Do ind=1, ind_max
                  write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
               End Do
               write(lun,*)

               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),1)
                  End Do
                  write(lun,*)
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),2)
                  End Do
                  write(lun,*)
               End Do
               write(lun,*)

            End Do
         End Do
         if ( iksnl /= (kv3/nspin) ) write(lun,*)
      End Do
      deallocate(norm2); deallocate(softpart)
      deallocate( compr_wk );      deallocate( compi_wk )
    end subroutine write_file_col_fbz

    subroutine write_file_col_fbz_sort( lun, is )
      integer, intent(in) :: lun, is

      integer :: ik, iksnl, ia, it, ib, jksnl, ii, iktmp
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind, lmtt2
      real(kind=DP) :: c1, weight, sign1
      complex(kind=CMPLDP) :: z1
      real(kind=DP), allocatable :: norm2(:,:), softpart(:,:,:)
      real(DP),allocatable :: crotylm(:,:,:)
      integer :: ikstart, ikend, jk, iopr, trev, mm, jorb
      real(kind=DP), allocatable :: compr_wk(:), compi_wk(:)

      if ( mype ==0 .and. is==1 ) then
         write(lun,'(A)') "PROJ  (orbital) lm decomposed + phase factor"
      endif

      allocate( norm2(natm+1,ind_max+1) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,2) );  softpart = 0.0d0
      if ( mype ==0 ) then
         write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3_fbz/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
         write(lun,*)
      endif

      if ( split_procar_file == ON ) then
         ikstart = ista_k +is -1;    ikend = iend_k
      else
         ikstart = is;               ikend = kv3
      endif

      allocate( compr_wk(nlmta_phi) )
      allocate( compi_wk(nlmta_phi) )

      Do jk=is, kv3_fbz, nspin

         iopr = iopr_k_fbz_to_ibz(jk)
         trev = trev_k_fbz_to_ibz(jk)
         jksnl = (jk-1)/nspin +1

         ikloop: Do ik=is, kv3, nspin
            Do ii=1, num_star_of_k(ik)
               if ( jk == star_of_k(ik,ii) ) then
                  iktmp = ik;  exit ikloop
               endif
            End Do
         end Do ikloop
         ik = iktmp

         iksnl = (ik-1)/nspin + 1

         weight = kv3 *qwgt(ik) /dble(ndim_spinor)
         if ( nspin == 1 ) weight = weight /2.0d0

         write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
              &     " k-point", jksnl, " :", vkxyz_fbz(jk,1:3,BUCS), &
              &     " weight = ", 1.0d0 / (kv3_fbz/nspin)
         write(lun,*)

         Do ib=1, neg
            softpart = 0.0d0;  norm2 = 0.0d0

            write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                 &     "band ", ib, &
                 &     " # energy ", ( eko_wk(ib,ik) -efermi )*Hartree, &
                 &     " # occ.", occ_wk(ib,ik)/weight
            write(lun,*)
               !
            if ( k_symmetry(ik) == GAMMA ) then
            else
               compr_wk = 0.d0;          compi_wk = 0.d0
               do iorb=1,nlmta_phi
                  do mm=1,nrorb(iorb,iopr)
                     jorb = irorb(mm,iorb,iopr)
                     compr_wk(iorb) = compr_wk(iorb) +compr(ib,jorb,1,ik) &
                          &          *crorb(mm,iorb,iopr)
                     compi_wk(iorb) = compi_wk(iorb) +compi(ib,jorb,1,ik) &
                          &          *crorb(mm,iorb,iopr)
                  end do
               end do
            endif
            !
            Do ia=1, natm
               it = ityp(ia)
               if (iproj_group(ia) == 0) cycle

               do lmt=1,ilmt_phi(it)
                  il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                  tau = taup_phi(lmt,it)
                  ind = (il-1)**2 +im
                  iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                  z1 = dcmplx( compr_wk(iorb), compi_wk(iorb) )

                  sign1 = 1.0d0
                  if ( trev == 1 ) sign1 = -sign1
#if 0
                  if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) sign1 = -sign1
#endif
                  if ( sign1 < .0d0 ) z1 = conjg(z1)
#if 1
                  if ( determinant_op(iopr) < 0 ) then
                     z1 = z1 *( (-1)**(il -1) )
                  endif
#endif

                  c1 = z1*conjg(z1)*( 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )

                  norm2(ia,ind) = c1
                  softpart(ia,ind,1) = real(z1)
                  softpart(ia,ind,2) = aimag(z1)
               end do
            End Do
!
            Do ia=1, natm
               c1 = 0.0d0
               Do ind=1, ind_max
                  c1 = c1 +norm2(ia,ind)
               End Do
               norm2(ia,ind_max+1) = c1
            End Do
!
            Do ind=1, ind_max+1
               c1 = 0.0d0
               Do ia=1, natm
                  c1 = c1 +norm2(ia,ind)
               End Do
               norm2(natm+1,ind) = c1
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*) "   tot"

            Do ia=1, natm
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind))
               End Do
               write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1)
               write(lun,*)
            End Do

            write(lun,'(A)',advance='no') "tot "
            Do ind=1, ind_max
               write(lun,'(F18.12)',advance='no') norm2(natm+1,index_p2v(ind))
            End Do
            write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1)
            write(lun,*)

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*)

            Do ia=1, natm
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),1)
               End Do
               write(lun,*)
               write(lun,'(I4)',advance='no') ia
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') softpart(ia,index_p2v(ind),2)
               End Do
               write(lun,*)
            End Do
            write(lun,*)
         End Do
         if ( jksnl /= (kv3_fbz/nspin) ) write(lun,*)
      End Do
      deallocate(norm2); deallocate(softpart)
      deallocate( compr_wk );      deallocate( compi_wk )
    end subroutine write_file_col_fbz_sort

    subroutine set_spinor_basis
      integer :: ia
      real(kind=DP) :: theta, phi, spn_quant_dir(3)
      real(kind=DP) :: mnorm, mx, my, mz

      if ( sw_fix_global_quantz_axis == ON ) then
         spn_quant_dir(1:3) = Global_Quantz_Axis_Fixed(1:3)
         theta = acos( spn_quant_dir(3) )
         phi = atan2( spn_quant_dir(2), spn_quant_dir(1) )

         Do ia=1, natm
            rotated_spinor(ia,1,1) =  exp( -zi *phi /2.0d0 ) *cos( theta /2.0d0 )
            rotated_spinor(ia,2,1) =  exp(  zi *phi /2.0d0 ) *sin( theta /2.0d0 )

            rotated_spinor(ia,1,2) = -exp( -zi *phi /2.0d0 ) *sin( theta /2.0d0 )
            rotated_spinor(ia,2,2) =  exp(  zi *phi /2.0d0 ) *cos( theta /2.0d0 )
         End Do
      else
         Do ia=1, natm
            mx = RhoMag_on_atom(ia,2)
            my = RhoMag_on_atom(ia,3)
            mz = RhoMag_on_atom(ia,4)
            mnorm = sqrt( mx**2 +my**2 +mz**2 )
            if ( mnorm < 1.0D-8 ) then
               spn_quant_dir(1) = 0.0d0
               spn_quant_dir(2) = 0.0d0
               spn_quant_dir(3) = 1.0d0

               theta = acos( spn_quant_dir(3) )
               phi = atan2( spn_quant_dir(2), spn_quant_dir(1) )
            else
               spn_quant_dir(1) = mx /mnorm
               spn_quant_dir(2) = my /mnorm
               spn_quant_dir(3) = mz /mnorm

               theta = acos( spn_quant_dir(3) )
               phi = atan2( spn_quant_dir(2), spn_quant_dir(1) )

               if ( spn_quant_dir(1)**2 +spn_quant_dir(2) < 1.0D-14 ) then
                  phi = 0.0d0
               endif
            endif

            rotated_spinor(ia,1,1) =  exp( -zi *phi /2.0d0 ) *cos( theta /2.0d0 )
            rotated_spinor(ia,2,1) =  exp(  zi *phi /2.0d0 ) *sin( theta /2.0d0 )

            rotated_spinor(ia,1,2) = -exp( -zi *phi /2.0d0 ) *sin( theta /2.0d0 )
            rotated_spinor(ia,2,2) =  exp(  zi *phi /2.0d0 ) *cos( theta /2.0d0 )
         End Do
      endif

    end subroutine set_spinor_basis

    subroutine write_file_noncl( lun )
      integer, intent(in) :: lun

      integer :: is, ik, iksnl, ia, it, ib
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind, m1, m2
      integer :: istmp, ikstart, ikend
      real(kind=DP) :: c1, c2, c3, c4, weight, vec1(3), vec2(3)
      complex(kind=CMPLDP) :: z1, z2, ztmp1, ztmp2, ztmp3, ztmp4
      complex(kind=CMPLDP) :: z1_u, z1_d, z2_u, z2_d

      real(kind=DP), allocatable :: norm2(:,:,:), softpart(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwork(:,:,:)

      if ( mype == 0 ) then
         if ( ekmode == OFF .or. nk_in_the_process == 1 ) then
            if ( procar_norm_zaxis_to_magdir == ON ) then
               write(lun,'(A)') &
                    &  "PROJ  (orbital) lm decomposed + phase factor (noncol,magdir//Z)"
            else
               write(lun,'(A)') &
                    &  "PROJ  (orbital) lm decomposed + phase factor (noncol)"
            endif
         endif
      endif

      allocate( norm2(natm+1,ind_max+1,ndim_magmom) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,ndim_spinor,2) );  softpart = 0.0d0

      if ( mype == 0 ) then
         if ( ekmode == OFF ) then
            write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
            write(lun,*)
         else
            if ( nk_in_the_process == 1 ) then
               write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3_ek/nspin,  &
                    &                              "# of bands: ", neg,  &
                    &                              "# of ions:", natm
               write(lun,*)
            endif
         endif
      endif

      allocate( zwork(natm,ind_max,ndim_spinor) )

      if ( split_procar_file == ON ) then
         ikstart = ista_k;    ikend = iend_k
      else
         ikstart = 1;         ikend = kv3
      endif

      Do ik=ikstart, ikend, ndim_spinor
         iksnl = (ik-1)/ndim_spinor + 1
         write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
              &     " k-point", iksnl, " :", vkxyz(ik,1:3,BUCS), &
              &     " weight = ", qwgt(ik)
         write(lun,*)

         weight = kv3 *qwgt(ik) /dble(ndim_spinor)

         Do ib=1, neg
            softpart = 0.0d0;  norm2 = 0.0d0

            write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                 &     "band ", ib, &
                 &     " # energy ", ( eko_wk(ib,iksnl) -efermi )*Hartree, &
                 &     " # occ.", occ_wk(ib,iksnl)/weight
            write(lun,*)
            Do ia=1, natm
               it = ityp(ia)
               if(iproj_group(ia) == 0) cycle

               do lmt=1,ilmt_phi(it)
                  il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                  tau = taup_phi(lmt,it)

                  ind = (il-1)**2 +im
                  iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                  z1 = dcmplx( compr(ib,iorb,1,ik),   compi(ib,iorb,1,ik) )
                  z2 = dcmplx( compr(ib,iorb,1,ik+1), compi(ib,iorb,1,ik+1) )

                  ztmp1 =   z1*conjg(z1) +z2*conjg(z2)
                  ztmp2 =   z1*conjg(z2) +z2*conjg(z1)
                  ztmp3 = (-z1*conjg(z2) +z2*conjg(z1) )*zi
                  ztmp4 =   z1*conjg(z1) -z2*conjg(z2)

                  c1 = 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl)
                  norm2(ia,ind,1) = ztmp1 *c1
                  norm2(ia,ind,2) = ztmp2 *c1
                  norm2(ia,ind,3) = ztmp3 *c1
                  norm2(ia,ind,4) = ztmp4 *c1

                  softpart(ia,ind,1,1) = compr(ib,iorb,1,ik)
                  softpart(ia,ind,1,2) = compi(ib,iorb,1,ik)
                  softpart(ia,ind,2,1) = compr(ib,iorb,1,ik+1)
                  softpart(ia,ind,2,2) = compi(ib,iorb,1,ik+1)
               End Do
            End Do
!!
            if ( procar_norm_zaxis_to_magdir == ON ) then
               Do ia=1, natm
                  il = 2         ! behave like p-orbital
                  Do ind=1, ind_max
                     vec1(1:3) = norm2(ia,ind,2:4)
                     Do m1=1, 3
                        c1 = 0.0d0
                        Do m2=1, 3
                           c1 = c1 +rotated_sph(ia,il,m1,m2) *vec1(m2)
                        End Do
                        vec2(m1) = c1
                     End Do
                     norm2(ia,ind,2:4) = vec2(1:3)
                  End Do
               End Do
            endif
!
            Do ia=1, natm
               c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
               Do ind=1, ind_max
                  c1 = c1 +norm2(ia,ind,1)
                  c2 = c2 +norm2(ia,ind,2)
                  c3 = c3 +norm2(ia,ind,3)
                  c4 = c4 +norm2(ia,ind,4)
               End Do
               norm2(ia,ind_max+1,1) = c1
               norm2(ia,ind_max+1,2) = c2
               norm2(ia,ind_max+1,3) = c3
               norm2(ia,ind_max+1,4) = c4
            End Do
!
            Do ind=1, ind_max+1
               c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
               Do ia=1, natm
                  c1 = c1 +norm2(ia,ind,1)
                  c2 = c2 +norm2(ia,ind,2)
                  c3 = c3 +norm2(ia,ind,3)
                  c4 = c4 +norm2(ia,ind,4)
               End Do
               norm2(natm+1,ind,1) = c1
               norm2(natm+1,ind,2) = c2
               norm2(natm+1,ind,3) = c3
               norm2(natm+1,ind,4) = c4
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*) "   tot"

            Do istmp=1, ndim_magmom
               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind),istmp)
                  End Do
                  write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1,istmp)
                  write(lun,*)
               End Do
               write(lun,'(A)',advance='no') "tot "
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') &
                       &         norm2(natm+1,index_p2v(ind),istmp)
               End Do
               write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1,istmp)
               write(lun,*)
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*)

            zwork = 0.0d0

            Do ia=1, natm
               Do ind=1, ind_max
                  z1_u = dcmplx( softpart(ia,index_p2v(ind),1,1), &
                       &       softpart(ia,index_p2v(ind),1,2) )
                  z1_d = dcmplx( softpart(ia,index_p2v(ind),2,1), &
                       &       softpart(ia,index_p2v(ind),2,2) )
! 1st spinor basis
#if 1
                  z2_u = rotated_spinor(ia,1,1)
                  z2_d = rotated_spinor(ia,2,1)
#else
                  z2_u = 1.0d0
                  z2_d = 0.0d0
#endif
                  zwork(ia,ind,1) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
! 2nd spinor basis
#if 1
                  z2_u = rotated_spinor(ia,1,2)
                  z2_d = rotated_spinor(ia,2,2)
#else
                  z2_u = 0.0d0
                  z2_d = 1.0d0
#endif
                  zwork(ia,ind,2) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
               End do
            End Do
! =====
            Do is=1, ndim_spinor
               if ( sw_fix_global_quantz_axis == ON ) then
                  if ( is==1 ) then
                     write(lun,'(A)') "(globally major) "
                  else
                     write(lun,'(A)') "(globally minor) "
                  endif
               else
                  if ( is==1 ) then
                     write(lun,'(A)') "(locally major) "
                  else
                     write(lun,'(A)') "(locally minor) "
                  endif
               endif

               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') &
                          &       real(zwork(ia,ind,is))
                  End Do
                  write(lun,*)
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') &
                          &       aimag(zwork(ia,ind,is))
                  End Do
                  write(lun,*)
               End Do
! =====
            End Do
            write(lun,*)
         End Do
         if ( iksnl /= (kv3/nspin) ) write(lun,*)
      End Do
      deallocate( zwork )

      deallocate(norm2); deallocate(softpart)
    end subroutine write_file_noncl

    subroutine get_op_spinor( nopr, op_spinor )
      integer, intent(in) :: nopr
      complex(kind=CMPLDP) :: op_spinor(2,2,nopr)

      character(len=9) :: system
      integer :: nsym, i, j, k, l
      real(kind=DP) :: csum
      real(kind=DP), allocatable :: rot_space(:,:,:)
      complex(kind=CMPLDP), allocatable :: rot_spinor(:,:,:)

      if ( lattice_system_from_m_CS_SG == "hexagonal" ) then
         system = "hexagonal";    nsym = 24;
      else
         system = "cubic";        nsym = 48;
      endif

      allocate( rot_space(3,3,nsym) )
      allocate( rot_spinor(2,2,nsym) )
      call m_CS_gen_opr_rspace_full( system, nsym, rot_space )
      call m_CS_gen_opr_spinor_full( system, nsym, rot_spinor )
!
      Do i=1, nopr
         jloop: Do j=1, nsym
            csum = 0.0d0
            Do k=1, 3
               Do l=1, 3
                  csum = csum +( op(k,l,i) -rot_space(k,l,j) )**2
               End Do
            End Do
            if ( csum < 1.0D-4 ) then
               op_spinor(:,:,i) = rot_spinor(:,:,j)
               exit jloop
            endif
         End Do jloop
      End Do
      deallocate( rot_space );    deallocate( rot_spinor )
    end subroutine get_op_spinor

    subroutine write_file_noncl_fbz( lun )
      integer, intent(in) :: lun

      integer :: ik, iksnl, ia, it, ib, jksnl, is1, is2
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind, ii, m1, m2
      integer :: istmp, ikstart, ikend, jk, iopr, trev, mm, jorb
      real(kind=DP) :: c1, c2, c3, c4, weight, sign1, vec1(3), vec2(3)
      complex(kind=CMPLDP) :: z1, z2, ztmp1, ztmp2, ztmp3, ztmp4
      complex(kind=CMPLDP) :: z1_u, z1_d, z2_u, z2_d
      real(kind=DP), allocatable :: norm2(:,:,:), softpart(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwork(:,:,:), op_spinor(:,:,:)
      complex(kind=CMPLDP), allocatable :: zcomp_wk(:,:)

      if ( mype == 0 ) then
         if ( procar_norm_zaxis_to_magdir == ON ) then
            write(lun,'(A)') &
                 &  "PROJ  (orbital) lm decomposed + phase factor (noncol,magdir//Z)"
         else
            write(lun,'(A)') &
                 &  "PROJ  (orbital) lm decomposed + phase factor (noncol)"
         endif
      endif

      allocate( norm2(natm+1,ind_max+1,ndim_magmom) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,ndim_spinor,2) );  softpart = 0.0d0

      if ( mype == 0 ) then
         write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3_fbz/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
         write(lun,*)
      endif

      allocate( zwork(natm,ind_max,ndim_spinor) )

      allocate( op_spinor(2,2,nopr) )
      call get_op_spinor( nopr, op_spinor )

      if ( split_procar_file == ON ) then
         ikstart = ista_k;    ikend = iend_k
      else
         ikstart = 1;         ikend = kv3
      endif

      allocate( zcomp_wk(nlmta_phi,ndim_spinor) )

      Do ik=ikstart, ikend, ndim_spinor
         iksnl = (ik-1)/ndim_spinor + 1
         weight = kv3 *qwgt(ik) /dble(ndim_spinor)

         Do ii=1, num_star_of_k(ik)
            jk = star_of_k(ik,ii)
            iopr = iopr_k_fbz_to_ibz(jk)
            trev = trev_k_fbz_to_ibz(jk)
            jksnl = (jk-1)/ndim_spinor + 1

            write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
              &     " k-point", jksnl, " :", vkxyz_fbz(jk,1:3,BUCS), &
              &     " weight = ", 1.0d0 / (kv3_fbz/ndim_spinor)
            write(lun,*)

            Do ib=1, neg
               softpart = 0.0d0;  norm2 = 0.0d0

               write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                    &     "band ", ib, &
                    &     " # energy ", ( eko_wk(ib,iksnl) -efermi )*Hartree, &
                    &     " # occ.", occ_wk(ib,iksnl)/weight
               write(lun,*)

               zcomp_wk = 0.d0;
               do is1=1, ndim_spinor
                  do iorb=1,nlmta_phi
                     Do is2=1, ndim_spinor
                        do mm=1,nrorb(iorb,iopr)
                           jorb = irorb(mm,iorb,iopr)
#if 1
                           z1 = dcmplx( compr(ib,jorb,1,ik+is2-1),&
                                &       compi(ib,jorb,1,ik+is2-1) )
                           zcomp_wk(iorb,is1) = zcomp_wk(iorb,is1) &
                                &             +z1 *crorb(mm,iorb,iopr) &
                                &             *op_spinor(is1,is2,iopr)
!                               &             *op_spinor(is2,is1,iopr)
#else
                           z1 = dcmplx( compr(ib,iorb,1,ik+is2-1),&
                                &       compi(ib,iorb,1,ik+is2-1) )
                           zcomp_wk(jorb,is1) = zcomp_wk(jorb,is1) &
                                &             +z1 *crorb(mm,iorb,iopr) &
                                &             *op_spinor(is1,is2,iopr)
#endif
                        end do
                     end do
                  end do
               end do
!!               if ( trev == 1 ) zcomp_wk = dconjg( zcomp_wk )   ! ????

               Do ia=1, natm
                  it = ityp(ia)
                  if(iproj_group(ia) == 0) cycle

                  do lmt=1,ilmt_phi(it)
                     il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                     tau = taup_phi(lmt,it)

                     ind = (il-1)**2 +im
                     iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                     sign1 = 1.0d0
                     if ( trev == 1 ) sign1 = -sign1
                     if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) sign1 = -sign1

                     if ( sign1 > 0.0d0 ) then
                        z1 = zcomp_wk(iorb,1)
                        z2 = zcomp_wk(iorb,2)
                     else
                        z1 = -conjg( zcomp_wk(iorb,2) )*zi
                        z2 =  conjg( zcomp_wk(iorb,1) )*zi
                     endif
#if 1
                     if ( determinant_op(iopr) < 0 ) then
                        z1 = z1 *( (-1)**(il -1) )
                        z2 = z2 *( (-1)**(il -1) )
                     endif
#endif
                     ztmp1 =   z1*conjg(z1) +z2*conjg(z2)
                     ztmp2 =   z1*conjg(z2) +z2*conjg(z1)
                     ztmp3 = (-z1*conjg(z2) +z2*conjg(z1) )*zi
                     ztmp4 =   z1*conjg(z1) -z2*conjg(z2)

                     c1 = 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl)
                     norm2(ia,ind,1) = ztmp1 *c1
                     norm2(ia,ind,2) = ztmp2 *c1
                     norm2(ia,ind,3) = ztmp3 *c1
                     norm2(ia,ind,4) = ztmp4 *c1

                     softpart(ia,ind,1,1) = real( z1 )
                     softpart(ia,ind,1,2) = aimag( z1 )
                     softpart(ia,ind,2,1) = real( z2 )
                     softpart(ia,ind,2,2) = aimag( z2 )
                  End Do
               End Do
!
               if ( procar_norm_zaxis_to_magdir == ON ) then
                  Do ia=1, natm
                     il = 2         ! behave like p-orbital
                     Do ind=1, ind_max
                        vec1(1:3) = norm2(ia,ind,2:4)
                        Do m1=1, 3
                           c1 = 0.0d0
                           Do m2=1, 3
                              c1 = c1 +rotated_sph(ia,il,m1,m2) *vec1(m2)
                           End Do
                           vec2(m1) = c1
                        End Do
                        norm2(ia,ind,2:4) = vec2(1:3)
                     End Do
                  End Do
               endif

               Do ia=1, natm
                  c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
                  Do ind=1, ind_max
                     c1 = c1 +norm2(ia,ind,1)
                     c2 = c2 +norm2(ia,ind,2)
                     c3 = c3 +norm2(ia,ind,3)
                     c4 = c4 +norm2(ia,ind,4)
                  End Do
                  norm2(ia,ind_max+1,1) = c1
                  norm2(ia,ind_max+1,2) = c2
                  norm2(ia,ind_max+1,3) = c3
                  norm2(ia,ind_max+1,4) = c4
               End Do
!
               Do ind=1, ind_max+1
                  c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
                  Do ia=1, natm
                     c1 = c1 +norm2(ia,ind,1)
                     c2 = c2 +norm2(ia,ind,2)
                     c3 = c3 +norm2(ia,ind,3)
                     c4 = c4 +norm2(ia,ind,4)
                  End Do
                  norm2(natm+1,ind,1) = c1
                  norm2(natm+1,ind,2) = c2
                  norm2(natm+1,ind,3) = c3
                  norm2(natm+1,ind,4) = c4
               End Do

               write(lun,'(A)',advance='no') "ion "
               Do ind=1, ind_max
                  write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
               End Do
               write(lun,*) "   tot"

               Do istmp=1, ndim_magmom
                  Do ia=1, natm
                     write(lun,'(I4)',advance='no') ia
                     Do ind=1, ind_max
                        write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind),istmp)
                     End Do
                     write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1,istmp)
                     write(lun,*)
                  End Do
                  write(lun,'(A)',advance='no') "tot "
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') &
                          &         norm2(natm+1,index_p2v(ind),istmp)
                  End Do
                  write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1,istmp)
                  write(lun,*)
               End Do

               write(lun,'(A)',advance='no') "ion "
               Do ind=1, ind_max
                  write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
               End Do
               write(lun,*)

               zwork = 0.0d0

               Do ia=1, natm
                  Do ind=1, ind_max
                     z1_u = dcmplx( softpart(ia,index_p2v(ind),1,1), &
                          &       softpart(ia,index_p2v(ind),1,2) )
                     z1_d = dcmplx( softpart(ia,index_p2v(ind),2,1), &
                          &       softpart(ia,index_p2v(ind),2,2) )
! 1st spinor basis
#if 1
                     z2_u = rotated_spinor(ia,1,1)
                     z2_d = rotated_spinor(ia,2,1)
#else
                     z2_u = 1.0d0
                     z2_d = 0.0d0
#endif
                     zwork(ia,ind,1) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
! 2nd spinor basis
#if 1
                     z2_u = rotated_spinor(ia,1,2)
                     z2_d = rotated_spinor(ia,2,2)
#else
                     z2_u = 0.0d0
                     z2_d = 1.0d0
#endif
                     zwork(ia,ind,2) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
                  End do
               End Do
! =====
               Do is1=1, ndim_spinor
                  if ( sw_fix_global_quantz_axis == ON ) then
                     if ( is1==1 ) then
                        write(lun,'(A)') "(globally major) "
                     else
                        write(lun,'(A)') "(globally minor) "
                     endif
                  else
                     if ( is1==1 ) then
                        write(lun,'(A)') "(locally major) "
                     else
                        write(lun,'(A)') "(locally minor) "
                     endif
                  endif

                  Do ia=1, natm
                     write(lun,'(I4)',advance='no') ia
                     Do ind=1, ind_max
                        write(lun,'(F18.12)',advance='no') &
                             &       real(zwork(ia,ind,is1))
                     End Do
                     write(lun,*)
                     write(lun,'(I4)',advance='no') ia
                     Do ind=1, ind_max
                        write(lun,'(F18.12)',advance='no') &
                             &       aimag(zwork(ia,ind,is1))
                     End Do
                     write(lun,*)
                  End Do
                  ! =====
               End Do
               write(lun,*)
            End Do
         End Do
         if ( iksnl /= (kv3/nspin) ) write(lun,*)
      End Do

      deallocate( zwork )
      deallocate( op_spinor )
      deallocate(norm2); deallocate(softpart)
    end subroutine write_file_noncl_fbz

    subroutine write_file_noncl_fbz_sort( lun )
      integer, intent(in) :: lun

      integer :: ik, iksnl, ia, it, ib, jksnl, is1, is2, ii
      integer :: lmt, il, im, tau, ip, iorb, lmtt, ind, iktmp
      integer :: istmp, ikstart, ikend, jk, iopr, trev, mm, jorb
      real(kind=DP) :: c1, c2, c3, c4, weight, sign1
      complex(kind=CMPLDP) :: z1, z2, ztmp1, ztmp2, ztmp3, ztmp4
      complex(kind=CMPLDP) :: z1_u, z1_d, z2_u, z2_d
      real(kind=DP), allocatable :: norm2(:,:,:), softpart(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwork(:,:,:), op_spinor(:,:,:)
      complex(kind=CMPLDP), allocatable :: zcomp_wk(:,:)
      integer :: ipriprocar = 1

      if ( mype == 0 ) then
         if ( procar_norm_zaxis_to_magdir == ON ) then
            write(lun,'(A)') &
                 &  "PROJ  (orbital) lm decomposed + phase factor (noncol,magdir//Z)"
         else
            write(lun,'(A)') &
                 &  "PROJ  (orbital) lm decomposed + phase factor (noncol)"
         endif
      endif

      allocate( norm2(natm+1,ind_max+1,ndim_magmom) );    norm2 = 0.0d0
      allocate( softpart(natm,ind_max,ndim_spinor,2) );  softpart = 0.0d0

      if ( mype == 0 ) then
         write(lun,'(A,I7,5X,A,I7,5X,A,I7)') "# of k-points: ", kv3_fbz/nspin,  &
              &                              "# of bands: ", neg,  &
              &                              "# of ions:", natm
         write(lun,*)
      endif

      allocate( zwork(natm,ind_max,ndim_spinor) )

      allocate( op_spinor(2,2,nopr) )
      call get_op_spinor( nopr, op_spinor )
!
      if ( split_procar_file == ON ) then
         ikstart = ista_k;    ikend = iend_k
      else
         ikstart = 1;         ikend = kv3
      endif

      allocate( zcomp_wk(nlmta_phi,ndim_spinor) )

      Do jk=1, kv3_fbz, ndim_spinor

         iopr = iopr_k_fbz_to_ibz(jk)
         trev = trev_k_fbz_to_ibz(jk)
         jksnl = (jk-1)/ndim_spinor + 1

         ikloop: Do ik=1, kv3, ndim_spinor
            Do ii=1, num_star_of_k(ik)
               if ( jk == star_of_k(ik,ii) ) then
                  iktmp = ik;  exit ikloop
               endif
            End Do
         end Do ikloop
         ik = iktmp

         iksnl = (ik-1)/ndim_spinor + 1
         weight = kv3 *qwgt(ik) /dble(ndim_spinor)

         if ( ipriprocar > 1 ) then
            if ( mype == 0 ) then
               write(nfout,*) "jksnl :", jksnl, iksnl, iopr, trev, &
                    &     magmom_dir_inversion_opr_flag(iopr), determinant_op(iopr)
            endif
         endif

         write(lun,'(A,I7,A,3F16.12,3X,A,F16.12)') &
              &     " k-point", jksnl, " :", vkxyz_fbz(jk,1:3,BUCS), &
              &     " weight = ", 1.0d0 / (kv3_fbz/ndim_spinor)
         write(lun,*)

         Do ib=1, neg
            softpart = 0.0d0;  norm2 = 0.0d0

            write(lun,'(A,I5,A,F16.12,A,F16.12)') &
                 &     "band ", ib, &
                 &     " # energy ", ( eko_wk(ib,iksnl) -efermi )*Hartree, &
                 &     " # occ.", occ_wk(ib,iksnl)/weight
            write(lun,*)

            zcomp_wk = 0.d0;
            do is1=1, ndim_spinor
               do iorb=1,nlmta_phi
                  Do is2=1, ndim_spinor
                     do mm=1,nrorb(iorb,iopr)
                        jorb = irorb(mm,iorb,iopr)
#if 1
                        z1 = dcmplx( compr(ib,jorb,1,ik+is2-1),&
                                &       compi(ib,jorb,1,ik+is2-1) )
                        zcomp_wk(iorb,is1) = zcomp_wk(iorb,is1) &
                             &             +z1 *crorb(mm,iorb,iopr) &
                             &             *op_spinor(is1,is2,iopr)
!                             &             *op_spinor(is2,is1,iopr)
#else
                        z1 = dcmplx( compr(ib,iorb,1,ik+is2-1),&
                             &       compi(ib,iorb,1,ik+is2-1) )
                        zcomp_wk(jorb,is1) = zcomp_wk(jorb,is1) &
                             &             +z1 *crorb(mm,iorb,iopr) &
                             &             *op_spinor(is1,is2,iopr)
#endif
                     end do
                  end do
               end do
            end do
!!            if ( trev == 1 ) zcomp_wk = dconjg( zcomp_wk )   ! ????

            Do ia=1, natm
               it = ityp(ia)
               if(iproj_group(ia) == 0) cycle

               do lmt=1,ilmt_phi(it)
                  il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
                  tau = taup_phi(lmt,it)

                  ind = (il-1)**2 +im
                  iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

                  sign1 = 1.0d0
                  if ( trev == 1 ) sign1 = -sign1
                  if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) sign1 = -sign1

                  if ( sign1 > 0.0d0 ) then
                     z1 = zcomp_wk(iorb,1)
                     z2 = zcomp_wk(iorb,2)
                  else
                     z1 = -conjg( zcomp_wk(iorb,2) )*zi
                     z2 =  conjg( zcomp_wk(iorb,1) )*zi
                  endif
#if 1
                  if ( determinant_op(iopr) < 0 ) then
                     z1 = z1 *( (-1)**(il -1) )
                     z2 = z2 *( (-1)**(il -1) )
                  endif
#endif
                  ztmp1 =   z1*conjg(z1) +z2*conjg(z2)
                  ztmp2 =   z1*conjg(z2) +z2*conjg(z1)
                  ztmp3 = (-z1*conjg(z2) +z2*conjg(z1) )*zi
                  ztmp4 =   z1*conjg(z1) -z2*conjg(z2)
!
                  c1 = 1.0d0 +qorb(iorb)/norm_phig_mpi(lmtt,iksnl)
                  norm2(ia,ind,1) = ztmp1 *c1
                  norm2(ia,ind,2) = ztmp2 *c1
                  norm2(ia,ind,3) = ztmp3 *c1
                  norm2(ia,ind,4) = ztmp4 *c1

                  softpart(ia,ind,1,1) = real( z1 )
                  softpart(ia,ind,1,2) = aimag( z1 )
                  softpart(ia,ind,2,1) = real( z2 )
                  softpart(ia,ind,2,2) = aimag( z2 )
               End Do
            End Do
            !
            Do ia=1, natm
               c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
               Do ind=1, ind_max
                  c1 = c1 +norm2(ia,ind,1)
                  c2 = c2 +norm2(ia,ind,2)
                  c3 = c3 +norm2(ia,ind,3)
                  c4 = c4 +norm2(ia,ind,4)
               End Do
               norm2(ia,ind_max+1,1) = c1
               norm2(ia,ind_max+1,2) = c2
               norm2(ia,ind_max+1,3) = c3
               norm2(ia,ind_max+1,4) = c4
            End Do
            !
            Do ind=1, ind_max+1
               c1 = 0.0d0;   c2 = 0.0d0;    c3 = 0.0d0;   c4 = 0.0d0
               Do ia=1, natm
                  c1 = c1 +norm2(ia,ind,1)
                  c2 = c2 +norm2(ia,ind,2)
                  c3 = c3 +norm2(ia,ind,3)
                  c4 = c4 +norm2(ia,ind,4)
               End Do
               norm2(natm+1,ind,1) = c1
               norm2(natm+1,ind,2) = c2
               norm2(natm+1,ind,3) = c3
               norm2(natm+1,ind,4) = c4
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*) "   tot"

            Do istmp=1, ndim_magmom
               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') norm2(ia,index_p2v(ind),istmp)
                  End Do
                  write(lun,'(F18.12)',advance='no') norm2(ia,ind_max+1,istmp)
                  write(lun,*)
               End Do
               write(lun,'(A)',advance='no') "tot "
               Do ind=1, ind_max
                  write(lun,'(F18.12)',advance='no') &
                       &         norm2(natm+1,index_p2v(ind),istmp)
               End Do
               write(lun,'(F18.12)',advance='no') norm2(natm+1,ind_max+1,istmp)
               write(lun,*)
            End Do

            write(lun,'(A)',advance='no') "ion "
            Do ind=1, ind_max
               write(lun,'(A,A)',advance='no') "    ", orb_name(index_p2v(ind))
            End Do
            write(lun,*)

            zwork = 0.0d0

            Do ia=1, natm
               Do ind=1, ind_max
                  z1_u = dcmplx( softpart(ia,index_p2v(ind),1,1), &
                       &       softpart(ia,index_p2v(ind),1,2) )
                  z1_d = dcmplx( softpart(ia,index_p2v(ind),2,1), &
                       &       softpart(ia,index_p2v(ind),2,2) )
                  ! 1st spinor basis
#if 1
                  z2_u = rotated_spinor(ia,1,1)
                  z2_d = rotated_spinor(ia,2,1)
#else
                  z2_u = 1.0d0
                  z2_d = 0.0d0
#endif
                  zwork(ia,ind,1) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
                  ! 2nd spinor basis
#if 1
                  z2_u = rotated_spinor(ia,1,2)
                  z2_d = rotated_spinor(ia,2,2)
#else
                  z2_u = 0.0d0
                  z2_d = 1.0d0
#endif
                  zwork(ia,ind,2) = conjg(z2_u) *z1_u +conjg(z2_d) *z1_d
               End do
            End Do
            ! =====
            Do is1=1, ndim_spinor
               if ( sw_fix_global_quantz_axis == ON ) then
                  if ( is1==1 ) then
                     write(lun,'(A)') "(globally major) "
                  else
                     write(lun,'(A)') "(globally minor) "
                  endif
               else
                  if ( is1==1 ) then
                     write(lun,'(A)') "(locally major) "
                  else
                     write(lun,'(A)') "(locally minor) "
                  endif
               endif

               Do ia=1, natm
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') &
                          &       real(zwork(ia,ind,is1))
                  End Do
                  write(lun,*)
                  write(lun,'(I4)',advance='no') ia
                  Do ind=1, ind_max
                     write(lun,'(F18.12)',advance='no') &
                          &       aimag(zwork(ia,ind,is1))
                  End Do
                  write(lun,*)
               End Do
               ! =====
            End Do
            write(lun,*)
         End Do
         if ( jksnl /= (kv3_fbz/nspin) ) write(lun,*)
      End Do
      deallocate( zwork )
      deallocate( op_spinor )
      deallocate(norm2); deallocate(softpart)
    end subroutine write_file_noncl_fbz_sort

    subroutine set_index_converter_p2v( ind_max, index_p2v )
      integer, intent(in) :: ind_max
      integer, intent(out) :: index_p2v(ind_max)
!
      index_p2v(1) = 1
      index_p2v(2) = 3
      index_p2v(3) = 4
      index_p2v(4) = 2
      index_p2v(5) = 7
      index_p2v(6) = 8
      index_p2v(7) = 5
      index_p2v(8) = 9
      index_p2v(9) = 6
    end subroutine set_index_converter_p2v

    subroutine set_orb_name( ind_max, orb_name )
      integer, intent(in) :: ind_max
      character*3, intent(out) :: orb_name(ind_max)

      if ( noncol .and. procar_sph_zaxis_to_magdir == ON ) then
         orb_name(1) = "  s";     orb_name(2) = " pX";       orb_name(3) = " pY"
         orb_name(4) = " pZ";     orb_name(5) = "dZ2";       orb_name(6) = "dX2"
         orb_name(7) = "dXY";     orb_name(8) = "dYZ";       orb_name(9) = "dZX"
      else
         orb_name(1) = "  s";     orb_name(2) = " px";       orb_name(3) = " py"
         orb_name(4) = " pz";     orb_name(5) = "dz2";       orb_name(6) = "dx2"
         orb_name(7) = "dxy";     orb_name(8) = "dyz";       orb_name(9) = "dzx"
      endif
    end subroutine set_orb_name

  end subroutine m_OP_wd_PROCAR

end module m_OP_decomp_band
