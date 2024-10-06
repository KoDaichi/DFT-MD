module m_OP_rotation
  use m_Control_Parameters,  only : nspin, ndim_chgpot, ndim_magmom, neg, ndim_spinor, &
       &                            noncol, YES, num_extra_bands, num_proj_elems, &
       &                            proj_group, proj_attribute, ipriorb_rot, &
       &                            sw_band_unfolding
  use m_Const_Parameters,    only : OFF, ON, DP, CMPLDP, GAMMA, BUCS, Bohr, PAI4, zi
  use m_Kpoints,     only : kv3, k_symmetry, kv3_ek, vkxyz, vkxyz_refcell
  use m_PseudoPotential,  only : nloc, ilmt_phi, lmta_phi, ltp_phi, mtp_phi, taup_phi, &
       &                         iproj_phi, lmtt_phi, qorb, nlmtt_phi, nlmta_phi, &
       &                          m_PP_tell_iorb_lmtt, m_PP_tell_iorb_ia_l_m_tau
  use m_Crystal_Structure,   only : nopr, altv
  use m_Ionic_System,    only : natm, ityp, iproj_group, if_pdos, speciesname, iatomn, &
       &                         cps, atom_key
  use m_NonLocal_Potential, only : norm_phig
  use m_Parallelization,    only : mype, map_k, map_e, myrank_k, myrank_e, map_z, npes, &
       &                           MPI_CommGroup, ierr, np_e, ista_k, iend_k, &
       &                           ista_e, iend_e, istep_e, nrank_e
  use m_Parallelization,    only : nrank_g
  use m_Electronic_Structure,  only : compr_l, compi_l, occup_l, neordr
  use m_IterationNumbers,      only : nk_in_the_process

  use m_ES_Noncollinear,   only : m_ES_set_Pauli_Matrix
  use m_SpinOrbit_Potential,  only : Mat_LS_with_real_ylm_L0, Mat_LS_with_real_ylm_L1, &
       &                             Mat_LS_with_real_ylm_L2, Mat_LS_with_real_ylm_L3, &
       &                             MatU_ylm_RC_L0, MatU_ylm_RC_L1, &
       &                             MatU_ylm_RC_L2, MatU_ylm_RC_L3
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP), allocatable :: porb_rot_matrix_real(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: porb_rot_matrix_cmplx(:,:,:,:)

  real(kind=DP), allocatable :: compr_rot_l(:,:,:,:), compi_rot_l(:,:,:,:)

  real(kind=DP), allocatable :: score_sigma_bond(:,:,:)
  character(len=8), allocatable :: orb_irrep(:,:,:)

  real(kind=DP), allocatable :: dm(:,:,:,:,:)
  complex(kind=CMPLDP), allocatable :: dm_ssrep(:,:,:,:,:)

  complex(kind=CMPLDP), allocatable :: Jcoeff_L0(:,:), JCoeff_L1(:,:)
  complex(kind=CMPLDP), allocatable :: Jcoeff_L2(:,:), JCoeff_L3(:,:)

contains

! ===========================================================-
  subroutine m_OP_calc_dm_col
    use m_Control_Parameters,  only : sw_write_orb_dens_mat_file
    use m_Files,  only :  m_Files_open_nfporb_dens_mat, &
         &                m_Files_close_nfporb_dens_mat, nfporb_dens_mat

    integer :: lmax, mmax

    lmax = nloc;   mmax = 2*lmax -1

    if ( allocated(dm) ) deallocate(dm)
    allocate( dm(natm,lmax,mmax,mmax,nspin) ); dm = 0.0d0
    call set_dm

    if ( sw_write_orb_dens_mat_file == ON ) call save_dm

  contains

    subroutine set_dm
      integer :: ik, is, iksnl, ia, it, iopr
      integer :: il, im, im1, im2, ip, i, i1, i2
      integer :: lmt, lmt1, lmt2, lmtt, lmtt1, lmtt2
      integer :: ib

      real(kind=DP) :: fac, c1, c2, temp
      real(kind=DP), allocatable :: dm_mpi(:,:,:,:,:)

      allocate( dm_mpi(natm,lmax,mmax,mmax,nspin) ); dm_mpi = 0.0d0

      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle

         iksnl = (ik-1)/nspin + 1
         is = mod(ik-1,nspin) +1

         do ib = 1, neg
            if(map_e(ib) == myrank_e) then
               do ia=1,natm
                  it = ityp(ia)
                  if(iproj_group(ia) == 0) cycle

                ! diagonal part
                  do lmt=1,ilmt_phi(it)
                     il = ltp_phi(lmt,it)
                     im  = mtp_phi(lmt,it)
                     ip  = iproj_phi(lmt,it)
                     i = lmta_phi(lmt,ia)
                     lmtt = lmtt_phi(lmt,it)

                     temp = 0.d0
                     if(k_symmetry(ik) == GAMMA) then
                        do iopr = 1,nopr
                           temp = temp + compr_l(map_z(ib),i,iopr,ik)**2 /2.0 &
                                &      *(1.d0+qorb(i)/(norm_phig(lmtt,iksnl)*2.0))
                        end do
                     else
                        do iopr = 1,nopr
                           temp = temp + (  compr_l(map_z(ib),i,iopr,ik)**2 &
                                &         + compi_l(map_z(ib),i,iopr,ik)**2) &
                                &      *(1.d0+qorb(i)/norm_phig(lmtt,iksnl))
                        end do
                     end if

                     dm_mpi(ia,il,im,im,is) = dm_mpi(ia,il,im,im,is) &
                          &                  + occup_l(map_z(ib),ik)*temp/dble(nopr)
                  end do

                  ! non-diagonal part
                  do lmt2=1,ilmt_phi(it)
                     do lmt1=1,lmt2-1
                        if(ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
                        il = ltp_phi(lmt1,it)
                        ip  = iproj_phi(lmt1,it)
                        im1  = mtp_phi(lmt1,it)
                        i1 = lmta_phi(lmt1,ia)
                        im2  = mtp_phi(lmt2,it)
                        i2 = lmta_phi(lmt2,ia)

                        lmtt1 = lmtt_phi(lmt1,it);     lmtt2 = lmtt_phi(lmt2,it)

                        c1 = ( qorb(i1) +qorb(i2) ) /2.0d0      ! approx
                        c2 = sqrt(norm_phig(lmtt1,iksnl)*norm_phig(lmtt2,iksnl))

                        temp = 0.d0

                        if(k_symmetry(ik) == GAMMA) then ! dame
                           do iopr = 1,nopr
                              temp = temp &
                                   & + compr_l(map_z(ib),i1,iopr,ik) &
                                   &  *compr_l(map_z(ib),i2,iopr,ik) /2.d0 &
                                   &  *( 1.0d0 +c1 /( c2 *2.d0) )
                           end do
                        else
                           do iopr = 1,nopr
                              temp = temp &
                                   & + ( compr_l(map_z(ib),i1,iopr,ik) &
                                   &    *compr_l(map_z(ib),i2,iopr,ik) &
                                   &    +compi_l(map_z(ib),i1,iopr,ik) &
                                   &    *compi_l(map_z(ib),i2,iopr,ik) )&
                                   &  *( 1.0d0 +c1 /c2 )
                           end do
                        end if

                        dm_mpi(ia,il,im1,im2,is) = dm_mpi(ia,il,im1,im2,is) &
                             &                   + occup_l(map_z(ib),ik)*temp/dble(nopr)
                        dm_mpi(ia,il,im2,im1,is) = dm_mpi(ia,il,im1,im2,is)
                     end do
                  end do

               end do
            end if
         end do
      end do

      if ( npes > 1 ) then
         call mpi_allreduce( dm_mpi, dm, natm *lmax *mmax**2 *nspin, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
      else
         dm = dm_mpi
      endif
      deallocate( dm_mpi )

      fac = 2.d0/kv3;
      dm = fac*dm
      if(nspin == 1) dm = 0.5d0*dm

    end subroutine set_dm

    subroutine save_dm
      integer :: ia, il, immax, im, im2, is
      integer :: lun

      if ( mype /= 0 ) return

      call m_Files_open_nfporb_dens_mat(0)

      lun = nfporb_dens_mat

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         write(lun) ia
         Do is=1, nspin
            write(lun) is
            Do il=1, nloc
               immax = 2*il -1
               write(lun) il
               Do im=1, immax
                  write(lun) ( dm(ia,il,im2,im,is), im2=1, immax )
               End Do
            End do
         End Do
      End Do
      call m_Files_close_nfporb_dens_mat

    end subroutine save_dm

  end subroutine m_OP_calc_dm_col

  subroutine m_OP_calc_dm_noncl         ! ssrep
    use m_Control_Parameters,  only : sw_write_orb_dens_mat_file
    use m_Files,  only :  m_Files_open_nfporb_dens_mat, &
         &                m_Files_close_nfporb_dens_mat, nfporb_dens_mat

    integer :: lmax, mmax

    lmax = nloc;   mmax = 2*lmax -1

    if ( allocated(dm_ssrep) ) deallocate(dm_ssrep)
    allocate( dm_ssrep(natm,lmax,mmax,mmax,ndim_chgpot) ); dm_ssrep = 0.0d0

    call set_dm_ssrep
    if ( sw_write_orb_dens_mat_file == ON ) call save_dm_ssrep

  contains

    subroutine set_dm_ssrep
      integer :: ik, is1, is2, iksnl, ia, it, iopr
      integer :: il, im, im1, im2, ip, i, i1, i2
      integer :: lmt, lmt1, lmt2, lmtt, lmtt1, lmtt2
      integer :: ib, istmp

      real(kind=DP) :: fac, c1, c2
      complex(kind=CMPLDP) :: ztemp, z1, z2
      complex(kind=CMPLDP), allocatable :: dm_mpi(:,:,:,:,:)

      allocate( dm_mpi(natm,lmax,mmax,mmax,ndim_chgpot) ); dm_mpi = 0.0d0

      do ik = 1, kv3, ndim_spinor
         if(map_k(ik) /= myrank_k) cycle

         iksnl = (ik-1)/nspin + 1

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               istmp = ( is1 -1 )*ndim_spinor + is2

               do ib = 1, neg
                  if(map_e(ib) == myrank_e) then
                     do ia=1,natm
                        it = ityp(ia)
                        if(iproj_group(ia) == 0) cycle

                ! diagonal part
                        do lmt=1,ilmt_phi(it)
                           il = ltp_phi(lmt,it)
                           im  = mtp_phi(lmt,it)
                           ip  = iproj_phi(lmt,it)
                           i = lmta_phi(lmt,ia)
                           lmtt = lmtt_phi(lmt,it)

                           ztemp = 0.0d0
                           do iopr=1,nopr
                              z1 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is1-1 ), &
                                   &       compi_l(map_z(ib),i,iopr,ik+is1-1 ) )
                              z2 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is2-1 ), &
                                   &       compi_l(map_z(ib),i,iopr,ik+is2-1 ) )
                              ztemp = ztemp + z1 *conjg(z2) &
                                   &      *( 1.d0+qorb(i)/norm_phig(lmtt,iksnl) )
                           end do

                           dm_mpi(ia,il,im,im,istmp) &
                                &    = dm_mpi(ia,il,im,im,istmp) &
                                &     + occup_l(map_z(ib),ik)*ztemp/dble(nopr)
                        end do

                  ! non-diagonal part
                        do lmt2=1,ilmt_phi(it)
!                           do lmt1=1,lmt2-1
                           do lmt1=1,ilmt_phi(it)
                              if ( lmt1 == lmt2 ) cycle
                              if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle

                              il = ltp_phi(lmt1,it)
                              ip  = iproj_phi(lmt1,it)
                              im1  = mtp_phi(lmt1,it)
                              i1 = lmta_phi(lmt1,ia)
                              im2  = mtp_phi(lmt2,it)
                              i2 = lmta_phi(lmt2,ia)

                              lmtt1 = lmtt_phi(lmt1,it);     lmtt2 = lmtt_phi(lmt2,it)

                              c1 = ( qorb(i1) +qorb(i2) ) /2.0d0      ! approx
                              c2 = sqrt(norm_phig(lmtt1,iksnl)*norm_phig(lmtt2,iksnl))

                              ztemp = 0.d0

                              do iopr = 1,nopr
                                 z1 = dcmplx( compr_l(map_z(ib),i1,iopr,ik+is1-1 ), &
                                      &       compi_l(map_z(ib),i1,iopr,ik+is1-1 ) )
                                 z2 = dcmplx( compr_l(map_z(ib),i2,iopr,ik+is2-1 ), &
                                      &       compi_l(map_z(ib),i2,iopr,ik+is2-1 ) )
                                 ztemp = ztemp + z1 *conjg(z2) &
                                      &  *( 1.0d0 +c1 /c2 )
                              end do

                              dm_mpi(ia,il,im1,im2,istmp) &
                                   &   = dm_mpi(ia,il,im1,im2,istmp) &
                                   &    + occup_l(map_z(ib),ik)*ztemp/dble(nopr)
!                              dm_mpi(ia,il,im2,im1,istmp) &
!                                   &  = conjg( dm_mpi(ia,il,im1,im2,istmp) )
                           end do
                        end do
                     end do
                  end if
               end do
            end Do
         end Do
      end do

      if ( npes > 1 ) then
         call mpi_allreduce( dm_mpi, dm_ssrep, natm *lmax *mmax**2 *ndim_chgpot *2, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
      else
         dm_ssrep = dm_mpi
      endif
      deallocate( dm_mpi )

      fac = 2.d0/kv3
      dm_ssrep = fac*dm_ssrep

    end subroutine set_dm_ssrep

    subroutine save_dm_ssrep
      integer :: ia, il, immax, im, im2, is
      integer :: lun

      if ( mype /= 0 ) return

      call m_Files_open_nfporb_dens_mat(0)

      lun = nfporb_dens_mat

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         write(lun) ia
         Do is=1, ndim_chgpot
            write(lun) is
            Do il=1, nloc
               immax = 2*il -1
               write(lun) il
               Do im=1, immax
                  write(lun) ( dm_ssrep(ia,il,im2,im,is), im2=1, immax )
               End Do
            End do
         End Do
      End Do
      call m_Files_close_nfporb_dens_mat

    end subroutine save_dm_ssrep

  end subroutine m_OP_calc_dm_noncl
  ! ===========================================================-

  subroutine m_OP_rotation_col(nfout)
    use m_Control_Parameters,  only : population_diag_mode, &
         &                            sw_calc_score_sigma_bond
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX, &
         &                            LOCAL_POINT_GROUP
    use m_Files,               only : m_Files_open_nfporb_rot_mat, &
         &                            m_Files_close_nfporb_rot_mat, nfporb_rot_mat
    integer, intent(in) :: nfout

    integer:: lmax, mmax
    real(kind=DP), allocatable :: porb2(:,:,:,:)

    lmax = nloc;   mmax = 2*lmax -1

    if ( allocated(porb_rot_matrix_real) ) deallocate( porb_rot_matrix_real )

    allocate( porb_rot_matrix_real(natm,lmax,mmax,mmax) );
    porb_rot_matrix_real = 0.0d0

    select case (population_diag_mode)
    case (DIAG_CHARGE_DENSITY_MATRIX)
       call diagonalization
    case (LOCAL_POINT_GROUP)
       allocate( orb_irrep(natm,lmax,mmax) )
       call m_OP_calc_symm
    end select

    if ( sw_calc_score_sigma_bond == ON ) then
       allocate( score_sigma_bond(natm,lmax,mmax) );   score_sigma_bond = -99.0d0
       call m_OP_calc_score_sigma_bond
    endif

    allocate( porb2(natm,lmax,mmax,nspin) ); porb2 = 0.0d0
    call set_porb2
    call printout
    call save_porb_rot_matrix_real

!    deallocate(dm);
    deallocate(porb2)
!    deallocate(score_sigma_bond)
!    if (allocated(orb_irrep)) deallocate(orb_irrep)

  contains

    subroutine save_porb_rot_matrix_real
      integer :: ia, il, immax, im, im2
      integer :: lun

      if ( mype /= 0 ) return

      call m_Files_open_nfporb_rot_mat(0)

      lun = nfporb_rot_mat
      write(lun) DIAG_CHARGE_DENSITY_MATRIX

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         write(lun) ia
         Do il=1, nloc
            immax = 2*il -1
            write(lun) il
            Do im=1, immax
               write(lun) ( porb_rot_matrix_real(ia,il,im2,im), im2=1, immax )
            End do
         End Do
      End Do
      call m_Files_close_nfporb_rot_mat

    end subroutine save_porb_rot_matrix_real

    subroutine printout
      integer :: ia, it, lmt, ilmta, il, im, tau, immax
      integer :: im2
      real(kind=DP) :: pup, pdn

      write(nfout,*)
      write(nfout,'(" --------- Orbital population rotated ---------")')

      if ( ipriorb_rot > 1 ) then
         write(nfout,'("  ia   l   m''  element coeff (1:2*l+1) ")')
         do ia=1,natm
            if (if_pdos(ia) == OFF) cycle
            if (iproj_group(ia) == 0) cycle
            it = ityp(ia)

            do lmt=1,ilmt_phi(it)
               ilmta = lmta_phi(lmt,ia)
               il = ltp_phi(lmt,it);   im = mtp_phi(lmt,it);   tau = taup_phi(lmt,it)
               immax = 2*il -1

               if ( immax <= 5 ) then
                  write(nfout,'(3(1x,i3),4x,a)',advance="no") &
                       &             ia, il-1, im, speciesname(it)
               else
                  write(nfout,'(3(1x,i3),4x,a)') ia, il-1, im, speciesname(it)
               endif
               Do im2=1, immax
                  write(nfout,'(1x,f10.5)',advance='no') &
                       &     ( porb_rot_matrix_real(ia,il,im2,im) )
               End Do
               write(nfout,*)
            end do
         end do
         write(nfout,*)
      endif

      pup=0.d0; pdn=0.d0

      write(nfout,'(A)',advance='no') &
           &    "  ia   l   m'  t    Porb(UP)   Porb(DN)  element"
      if ( sw_calc_score_sigma_bond == ON ) then
         write(nfout,'(A)',advance='no') "  score_bond"
      endif
      if ( population_diag_mode == LOCAL_POINT_GROUP ) then
         write(nfout,'(A)') "  symmetry"
      else
         write(nfout,*)
      endif

      do ia=1,natm
         if (if_pdos(ia) == OFF) cycle
         if (iproj_group(ia) == 0) cycle
         it = ityp(ia)

         do lmt=1,ilmt_phi(it)
            ilmta = lmta_phi(lmt,ia)
            il = ltp_phi(lmt,it);  im = mtp_phi(lmt,it);  tau = taup_phi(lmt,it)
            if ( nspin == 1 ) then
               write(nfout,'(4(1x,i3),2(1x,f10.5),4x,a4)',advance='no') &
                    &        ia, il-1, im, tau, &
                    &        porb2(ia,il,im,1), porb2(ia,il,im,1), &
                    &        speciesname(it)
               if ( sw_calc_score_sigma_bond == ON ) then
                  write(nfout,'(3x,f10.5)',advance='no') score_sigma_bond(ia,il,im)
               endif
               if ( population_diag_mode == LOCAL_POINT_GROUP ) then
                  write(nfout,'(4x,a8)') orb_irrep(ia,il,im)
               else
                  write(nfout,*)
               endif
               pup = pup +porb2(ia,il,im,1)

            else
               write(nfout,'(4(1x,i3),2(1x,f10.5),4x,a4)',advance='no') &
                    &        ia, il-1, im, tau, porb2(ia,il,im,1:2), &
                    &        speciesname(it)
               if ( sw_calc_score_sigma_bond == ON ) then
                  write(nfout,'(3x,f10.5)',advance='no') score_sigma_bond(ia,il,im)
               endif
               if ( population_diag_mode == LOCAL_POINT_GROUP ) then
                  write(nfout,'(4x,a8)') orb_irrep(ia,il,im)
               else
                  write(nfout,*)
               endif
               pup = pup +porb2(ia,il,im,1)
               pdn = pdn +porb2(ia,il,im,2)
            endif
         end do
      end do

      if ( nspin == 1 ) then
         pdn = pup
         write(nfout,'(8x,"Total : ",3(1x,f10.5))') pup,pdn,pup+pdn
      else
         write(nfout,'(8x,"Total : ",3(1x,f10.5))') pup,pdn,pup+pdn
      endif
      write(nfout,*)

    end subroutine printout

    subroutine set_porb2
      integer :: ia, il, immax, size2, im1, im2, im3, is
      real(kind=DP) :: ctmp, cwk(ndim_chgpot)

      do ia=1,natm
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1
            size2 = immax

            Do im1=1, size2
               DO is=1, ndim_magmom
                  ctmp = 0.0d0
                  Do im2=1, immax
                     Do im3=1, immax
                        ctmp = ctmp &
                             &    + porb_rot_matrix_real( ia,il,im2,im1 )  &
                             &     *dm(ia,il,im2,im3,is) &
                             &     *porb_rot_matrix_real( ia,il,im3,im1 )
                     End Do
                  End Do
                  cwk(is) = ctmp
               End DO
               porb2( ia, il, im1, 1:ndim_magmom ) = cwk(1:ndim_magmom)
            End Do
         End Do
      ENd do

    end subroutine set_porb2

    subroutine diagonalization
      integer :: is,ia,it,ig,i,ip,il, tau, im1, im2, im3
      integer :: lwork,info, immax
      real(kind=DP) :: abstol,nfound, ctmp

      integer, allocatable :: iwork(:),ifail(:)
      real(kind=DP), allocatable :: rwork(:), amat(:,:), evec(:), rvmat(:,:)
      real(kind=DP), external :: dlamch


      abstol = 2*dlamch('S')

     do ia=1,natm
         it = ityp(ia)
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1;
            lwork = 8 *immax
            allocate( rwork(lwork) );
            allocate( iwork(5*immax) );   allocate( ifail(immax) )
            allocate( amat(immax,immax) )
            allocate( evec(immax) )
            allocate( rvmat(immax,immax) )

            amat = 0.0d0
            Do im1=1, immax
               Do im2=1, immax
                  Do is=1, ndim_magmom
                     amat(im1,im2) = amat(im1,im2) +dm(ia,il,im1,im2,is)
                  End Do
               End Do
            ENd Do
            call dsyevx( 'V', 'A', 'U', immax, amat, immax, 0.d0, 0.d0, 0, 0, &
                 &       abstol, nfound, evec, rvmat, immax, rwork, lwork, &
                 &       iwork, ifail, info )
            if (info /= 0) then
               write(nfout,*) 'dsyevx: info=',info
            end if
            Do im1=1, immax
               Do im2=1, immax
                  porb_rot_matrix_real( ia,il,im2,im1 ) = rvmat(im2,im1)
               ENd Do
            End Do

            deallocate(rwork);  deallocate(iwork);   deallocate(ifail)
            deallocate(amat);   deallocate(evec);    deallocate(rvmat)
         end do
      end do
    end subroutine diagonalization
  end subroutine m_OP_rotation_col

  subroutine m_OP_rotation_noncl(nfout)
    use m_Control_Parameters,  only : population_diag_mode, &
         &                            sw_calc_score_sigma_bond
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX, &
         &                            DIAG_SPIN_DENSITY_MATRIX, &
         &                            DIAG_LS_with_t2g_octa, DIAG_LS, &
         &                            LOCAL_POINT_GROUP, LOCAL_DOUBLE_POINT_GROUP
    use m_Files,               only : m_Files_open_nfporb_rot_mat, &
         &                            m_Files_close_nfporb_rot_mat, nfporb_rot_mat

    integer, intent(in) :: nfout

    integer:: lmax, mmax

    real(kind=DP), allocatable :: porb2(:,:,:,:)
    real(kind=DP), allocatable :: dm(:,:,:,:,:)
!
    integer :: fac, i, j
!
    lmax = nloc;   mmax = 2*lmax -1

    if ( allocated(porb_rot_matrix_real) ) deallocate( porb_rot_matrix_real )
    if ( allocated(porb_rot_matrix_cmplx) ) deallocate( porb_rot_matrix_cmplx )

    select case (population_diag_mode)
    case (LOCAL_POINT_GROUP)
       allocate( dm(natm,lmax,mmax,mmax,ndim_magmom) ); dm = 0.0d0
       call set_dm

       allocate( porb_rot_matrix_real(natm,lmax,mmax,mmax) );
       porb_rot_matrix_real = 0.0d0

       allocate( orb_irrep(natm,lmax,mmax) )
       call m_OP_calc_symm
       if ( sw_calc_score_sigma_bond == ON ) then
          allocate( score_sigma_bond(natm,lmax,mmax) );   score_sigma_bond = -99.0d0
          call m_OP_calc_score_sigma_bond
       endif

       allocate( porb2(natm,lmax,mmax,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_chgmag_densmat
       call printout_chg_densmat_mode
       call save_porb_rot_matrix_real

    case (LOCAL_DOUBLE_POINT_GROUP)
       fac = ndim_spinor
       allocate( porb_rot_matrix_cmplx(natm,lmax,mmax*fac,mmax*fac) );
       porb_rot_matrix_cmplx = 0.0d0

       allocate( orb_irrep(natm,lmax,mmax*fac) )
       call m_OP_calc_symm

       allocate( porb2(natm,lmax,mmax*fac,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_spin_densmat
       call printout_spin_densmat_mode
       call save_porb_rot_matrix_cmplx

    case (DIAG_CHARGE_DENSITY_MATRIX)
       allocate( dm(natm,lmax,mmax,mmax,ndim_magmom) ); dm = 0.0d0
       call set_dm

       allocate( porb_rot_matrix_real(natm,lmax,mmax,mmax) );
       porb_rot_matrix_real = 0.0d0
       call diagonalize_chg_densmat
       if ( sw_calc_score_sigma_bond == ON ) then
          allocate( score_sigma_bond(natm,lmax,mmax) );   score_sigma_bond = -99.0d0
          call m_OP_calc_score_sigma_bond
       endif

       allocate( porb2(natm,lmax,mmax,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_chgmag_densmat
       call printout_chg_densmat_mode
       call save_porb_rot_matrix_real

    case (DIAG_SPIN_DENSITY_MATRIX)
       fac = ndim_spinor

       allocate( porb_rot_matrix_cmplx(natm,lmax,mmax*fac,mmax*fac) );
       porb_rot_matrix_cmplx = 0.0d0
       call diagonalize_spin_densmat
       if ( sw_calc_score_sigma_bond == ON ) then
          allocate( score_sigma_bond(natm,lmax,mmax*fac) );  score_sigma_bond = -99.0d0
          call m_OP_calc_score_sigma_bond
       endif

       allocate( porb2(natm,lmax,mmax*fac,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_spin_densmat
       call printout_spin_densmat_mode
       call save_porb_rot_matrix_cmplx

    case (DIAG_LS_with_t2g_OCTA)
       fac = ndim_spinor

       allocate( dm(natm,lmax,mmax,mmax,ndim_magmom) ); dm = 0.0d0
       call set_dm

       allocate( porb_rot_matrix_real(natm,lmax,mmax,mmax) );
       porb_rot_matrix_real = 0.0d0
       call diagonalize_chg_densmat

       allocate( porb_rot_matrix_cmplx(natm,lmax,mmax*fac,mmax*fac) );
       porb_rot_matrix_cmplx = 0.0d0
       call set_porb_rot_matrix_jeff

       allocate( porb2(natm,lmax,mmax*fac,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_spin_densmat
       call printout_jeff_mode
       call save_porb_rot_matrix_cmplx

    case (DIAG_LS)
       fac = ndim_spinor

       allocate( porb_rot_matrix_cmplx(natm,lmax,mmax*fac,mmax*fac) );
       porb_rot_matrix_cmplx = 0.0d0
       call set_porb_rot_matrix_jrep3

       allocate( porb2(natm,lmax,mmax*fac,ndim_magmom) ); porb2 = 0.0d0
       call set_porb2_from_spin_densmat
       call printout_jrep_mode3
       call save_porb_rot_matrix_cmplx

    end select

!    deallocate(dm_ssrep);
    deallocate(porb2)
    if ( allocated(dm) ) deallocate(dm)

  contains

    subroutine save_porb_rot_matrix_real
      integer :: ia, il, immax, im, im2
      integer :: lun

      if ( mype /= 0 ) return

      call m_Files_open_nfporb_rot_mat(0)

      lun = nfporb_rot_mat
      write(lun) population_diag_mode

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         write(lun) ia
         Do il=1, nloc
            immax = 2*il -1
            write(lun) il
            Do im=1, immax
               write(lun) ( porb_rot_matrix_real(ia,il,im2,im), im2=1, immax )
            End do
         End Do
      End Do
      call m_Files_close_nfporb_rot_mat
    end subroutine save_porb_rot_matrix_real

    subroutine save_porb_rot_matrix_cmplx
      integer :: ia, il, immax, im, im2
      integer :: lun

      if ( mype /= 0 ) return

      call m_Files_open_nfporb_rot_mat(0)

      lun = nfporb_rot_mat
      write(lun) population_diag_mode

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         write(lun) ia
         Do il=1, nloc
            immax = 2*il -1
            write(lun) il
            Do im=1, immax*2
               write(lun) ( porb_rot_matrix_cmplx(ia,il,im2,im), im2=1, immax*2 )
            End do
         End Do
      End Do
      call m_Files_close_nfporb_rot_mat
    end subroutine save_porb_rot_matrix_cmplx

    subroutine printout_chg_densmat_mode
      integer :: ia, it, lmt, ilmta, il, im, tau, immax
      integer :: im2

      write(nfout,*)
      write(nfout,'(" --------- Orbital population rotated ---------")')

      if ( ipriorb_rot > 1 ) then
         write(nfout,'("  ia   l   m''  element coeff (1:2*l+1) ")')
         do ia=1,natm
            if (if_pdos(ia) == OFF) cycle
            if (iproj_group(ia) == 0) cycle
            it = ityp(ia)

            do lmt=1,ilmt_phi(it)
               ilmta = lmta_phi(lmt,ia)
               il = ltp_phi(lmt,it);   im = mtp_phi(lmt,it);   tau = taup_phi(lmt,it)
               immax = 2*il -1

               write(nfout,'(3(1x,i3),4x,a)') ia, il-1, im, speciesname(it)
               write(nfout,'(7(1x,f10.5))') &
                    &     ( porb_rot_matrix_real(ia,il,im2,im), im2=1, immax )
               write(nfout,*)
            end do
         end do
      endif

      write(nfout,'(A)',advance='no') &
           &  "  ia   l   m'  t    Porb(tot)  Porb(mx)   Porb(my)   Porb(mz) element"
      if ( sw_calc_score_sigma_bond == ON ) then
         write(nfout,'(A)',advance='no') "  score"
      endif
      if ( population_diag_mode == LOCAL_POINT_GROUP ) then
         write(nfout,'(A)') "  symmetry"
      else
         write(nfout,*)
      endif

      do ia=1,natm
         if (if_pdos(ia) == OFF) cycle
         if (iproj_group(ia) == 0) cycle
         it = ityp(ia)

         do lmt=1,ilmt_phi(it)
            ilmta = lmta_phi(lmt,ia)
            il = ltp_phi(lmt,it);  im = mtp_phi(lmt,it);  tau = taup_phi(lmt,it)
            write(nfout,'(4(1x,i3),4(1x,f10.5),4x,a)',advance='no') &
                 &        ia, il-1, im, tau, porb2(ia,il,im,:),speciesname(it)
            if ( sw_calc_score_sigma_bond == ON ) then
               write(nfout,'(1x,f8.4)',advance='no') score_sigma_bond(ia,il,im)
            endif
            if ( population_diag_mode == LOCAL_POINT_GROUP ) then
               write(nfout,'(2x,a8)') orb_irrep(ia,il,im)
            else
               write(nfout,*)
            endif
         end do
      end do
      write(nfout,*)

    end subroutine printout_chg_densmat_mode

    subroutine printout_spin_densmat_mode
      integer :: ia, it, lmt, ilmta, il, im, tau, immax
      integer :: im1, im2, is1
      real(kind=DP) ::val_j

      write(nfout,*)
      write(nfout,'(" --------- Orbital population rotated ---------")')
      if ( ipriorb_rot > 1 ) then
         write(nfout,'(A)') "  ia   l  ms'   element coeff (1:(2*l+1)*2) "
         do ia=1,natm
            if (if_pdos(ia) == OFF) cycle
            if (iproj_group(ia) == 0) cycle
            it = ityp(ia)

            do lmt=1,ilmt_phi(it)
               ilmta = lmta_phi(lmt,ia)
               il = ltp_phi(lmt,it);   im = mtp_phi(lmt,it);   tau = taup_phi(lmt,it)
               immax = 2*il -1

               Do is1=1, ndim_spinor
                  im1 = ndim_spinor*(im-1) +is1

                  write(nfout,'(3(1x,i3),4x,a)') ia, il-1, im1, speciesname(it)
                  write(nfout,'(A)') "UP   spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Re: ", (real(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A)') "DOWN spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Re: ", (real(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                 im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                  im2=1, immax )
                  write(nfout,*)
               end Do
            end do
         end do
      endif

      write(nfout,'(A)',advance='no') &
           & "  ia   l  ms'  t    Porb(tot)  Porb(mx)   Porb(my)   Porb(mz) element"
      if ( population_diag_mode == LOCAL_DOUBLE_POINT_GROUP ) then
         write(nfout,'(A)',advance='no')  "  j"
         write(nfout,'(A)') "   symmetry"
      else
         if ( sw_calc_score_sigma_bond == ON ) then
            write(nfout,'(A)',advance='no') "  score"
         endif
         write(nfout,*)
      endif
      do ia=1,natm
         if (if_pdos(ia) == OFF) cycle
         if (iproj_group(ia) == 0) cycle
         it = ityp(ia)

         do lmt=1,ilmt_phi(it)
            ilmta = lmta_phi(lmt,ia)
            il = ltp_phi(lmt,it);  im = mtp_phi(lmt,it);  tau = taup_phi(lmt,it)
            Do is1=1, ndim_spinor
               im1 = ndim_spinor*(im-1) +is1
               write(nfout,'(4(1x,i3),4(1x,f10.5),4x,a)',advance='no') &
                    &        ia, il-1, im1, tau, porb2(ia,il,im1,:),speciesname(it)
               if ( population_diag_mode == LOCAL_DOUBLE_POINT_GROUP ) then
                  if ( im1 <= 2*(il-1) ) then
                     val_j = il -1 -0.5d0
                  else
                     val_j = il -1 +0.5d0
                  endif
                  write(nfout,'(2x,i1,a2)',advance='no') nint(2*val_j),  "/2"
                  write(nfout,'(2x,a8)') orb_irrep(ia,il,im1)
               else
                  if ( sw_calc_score_sigma_bond == ON ) then
                     write(nfout,'(1x,f8.4)',advance='no') score_sigma_bond(ia,il,im1)
                  endif
                  write(nfout,*)
               endif
            End Do
         end do
      end do
      write(nfout,*)

    end subroutine printout_spin_densmat_mode

    subroutine printout_jeff_mode
      integer :: ia, it, lmt, ilmta, il, im, tau, immax
      integer :: im1, im2, is1
      character*18 comment

      write(nfout,*)
      write(nfout,'(" --------- Orbital population rotated ---------")')

      write(nfout,'("  ia   l   ms''  element coeff (1:(2*l+1)*2) ")')
      if ( ipriorb_rot > 1 ) then
         do ia=1,natm
            if (if_pdos(ia) == OFF) cycle
            if (iproj_group(ia) == 0) cycle
            it = ityp(ia)

            do lmt=1,ilmt_phi(it)
               ilmta = lmta_phi(lmt,ia)
               il = ltp_phi(lmt,it);   im = mtp_phi(lmt,it);   tau = taup_phi(lmt,it)
               immax = 2*il -1
               if ( il /= 2+1 ) cycle

               Do is1=1, ndim_spinor
                  im1 = ndim_spinor*(im-1) +is1

                  comment = ""
                  if ( il == 2+1 ) then
                     if ( im1 <= 4 ) then
                        comment = "(t2g Jeff = 3/2)"
                     else if ( im1 <= 6 ) then
                        comment = "(t2g Jeff = 1/2)"
                     else
                        comment = "(eg)"
                     endif
                  endif

                  write(nfout,'(3(1x,i3),4x,a,12x,a)') ia, il-1, im1, &
                       &                               speciesname(it), comment
                  write(nfout,'(A)') "UP   spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Re: ", (real(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A)') "DOWN spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Re: ", (real(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                 im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                  im2=1, immax )
                  write(nfout,*)
               end Do
            end do
         end do
      endif

! ---
      write(nfout,'("  ia   l   ms''  t    Porb(tot)  Porb(mx)   Porb(my)   Porb(mz) element")')
      do ia=1,natm
         if (if_pdos(ia) == OFF) cycle
         if (iproj_group(ia) == 0) cycle
         it = ityp(ia)

         do lmt=1,ilmt_phi(it)
            ilmta = lmta_phi(lmt,ia)
            il = ltp_phi(lmt,it);  im = mtp_phi(lmt,it);  tau = taup_phi(lmt,it)
            if ( il /=2 +1 ) cycle

            Do is1=1, ndim_spinor
               im1 = ndim_spinor*(im-1) +is1

               comment = ""
               if ( il == 2+1 ) then
                  if ( im1 <= 4 ) then
                     comment = "(t2g Jeff = 3/2)"
                  else if ( im1 <= 6 ) then
                     comment = "(t2g Jeff = 1/2)"
                  else
                     comment = "(eg)"
                  endif
               endif
               write(nfout,'(4(1x,i3),4(1x,f10.5),4x,a,2x,a)') &
                    &        ia, il-1, im1, tau, porb2(ia,il,im1,:), &
                    &        speciesname(it), comment
            End Do
         end do
      end do
      write(nfout,*)

    end subroutine printout_jeff_mode

    subroutine printout_jrep_mode3
      integer :: ia, it, lmt, ilmta, il, im, tau, immax
      integer :: im1, im2, is1, my_l
      character*17 comment
      real(kind=DP) :: val_j, val_mj

      write(nfout,*)
      write(nfout,'(" --------- Orbital population rotated ---------")')

      if ( ipriorb_rot > 1 ) then
         write(nfout,'(A)') "  ia   l  ms'  element coeff (1:(2*l+1)*2) "
         do ia=1,natm
            if (if_pdos(ia) == OFF) cycle
            if (iproj_group(ia) == 0) cycle
            it = ityp(ia)

            do lmt=1,ilmt_phi(it)
               ilmta = lmta_phi(lmt,ia)
               il = ltp_phi(lmt,it);   im = mtp_phi(lmt,it);   tau = taup_phi(lmt,it)
               immax = 2*il -1

               Do is1=1, ndim_spinor
                  im1 = ndim_spinor*(im-1) +is1

                  my_l = il -1
                  if ( my_l > 0 ) then
                     if ( im1 <= 2*my_l ) then
                        val_j = my_l -0.5d0;     val_mj = -val_j +(im1 -1)
                     else
                        val_j = my_l +0.5d0;     val_mj = -val_j +(im1 -2*my_l -1)
                     endif
                  else
                     val_j = my_l +0.5d0;        val_mj = -val_j +(im1 -2*my_l -1)
                  endif
                  write(comment,'(A,I1,A,I2,A)') "(j=",nint(2*val_j),"/2, mj=", &
                       &                          nint(2*val_mj), "/2)"
                  write(nfout,'(3(1x,i3),4x,a,12x,a)') ia, il-1, im1, &
                       &                               speciesname(it), comment
                  write(nfout,'(A)') "UP   spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Re: ", (real(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &     "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,im2,im1)), &
                       &                   im2=1, immax )
                  write(nfout,'(A)') "DOWN spin component"
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Re: ", (real(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                 im2=1, immax )
                  write(nfout,'(A,7(1x,f10.5))') &
                       &   "Im: ", (aimag(porb_rot_matrix_cmplx(ia,il,immax+im2,im1)), &
                       &                  im2=1, immax )
                  write(nfout,*)
               end Do
            end do
         end do
      endif
! ---
      write(nfout,'(A,4X,A)') &
           &  "  ia   l  ms'  t    Porb(tot)  Porb(mx)   Porb(my)   Porb(mz) element", &
           &  "j     mj"
      do ia=1,natm
         if (if_pdos(ia) == OFF) cycle
         if (iproj_group(ia) == 0) cycle
         it = ityp(ia)

         do lmt=1,ilmt_phi(it)
            ilmta = lmta_phi(lmt,ia)
            il = ltp_phi(lmt,it);  im = mtp_phi(lmt,it);  tau = taup_phi(lmt,it)

            Do is1=1, ndim_spinor
               im1 = ndim_spinor*(im-1) +is1

               my_l = il -1
               if ( my_l > 0 ) then
                  if ( im1 <= 2*my_l ) then
                     val_j = my_l -0.5d0;       val_mj = -val_j +(im1 -1)
                  else
                     val_j = my_l +0.5d0;       val_mj = -val_j +(im1 -2*my_l -1)
                  endif
               else
                  val_j = my_l +0.5d0;          val_mj = -val_j +(im1 -2*my_l -1)
               endif
               write(comment,'(I1,A,I2,A)') nint(2*val_j),  "/2  ", &
                    &                       nint(2*val_mj), "/2"
               write(nfout,'(4(1x,i3),4(1x,f10.5),4x,a,2x,a,a)') &
                    &        ia, il-1, im1, tau, porb2(ia,il,im1,:), &
                    &        speciesname(it), "  ", adjustl(comment)
            End Do
         end do
      end do
      write(nfout,*)

    end subroutine printout_jrep_mode3

    subroutine set_porb_rot_matrix_jrep3
      integer :: ia, il, im1, im2, im3, im4, im0, is1, is2, istmp
      integer :: immax, lwork, info, size2
      complex(kind=CMPLDP) :: ztmp
      real(kind=DP), allocatable :: eigenval(:), rwork(:)
      complex(kind=CMPLDP), allocatable :: MatA(:,:), work(:), zwk(:,:)

      integer :: my_l, ma, m1, m2, num
      real(kind=DP) :: c1, c2, val_j, val_mj

      call set_Jcoeff_L0_L3

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         Do il=1, nloc
            immax = 2*il -1

            Do im0=1, immax *ndim_spinor
               Do im1=1, immax
                  Do is1=1, ndim_spinor
                     if ( il==1 ) then
                        porb_rot_matrix_cmplx(ia,il,immax*(is1-1)+im1,im0) &
                             & = JCoeff_L0(immax*(is1-1)+im1,im0)
                     else if ( il == 2 ) then
                        porb_rot_matrix_cmplx(ia,il,immax*(is1-1)+im1,im0) &
                             & = JCoeff_L1(immax*(is1-1)+im1,im0)
                     else if ( il == 3 ) then
                        porb_rot_matrix_cmplx(ia,il,immax*(is1-1)+im1,im0) &
                             & = JCoeff_L2(immax*(is1-1)+im1,im0)
                     else if ( il == 4 ) then
                        porb_rot_matrix_cmplx(ia,il,immax*(is1-1)+im1,im0) &
                             & = JCoeff_L3(immax*(is1-1)+im1,im0)
                     endif
                  End Do
               End Do
            End Do
         End Do
      End Do
      deallocate( JCoeff_L0 );      deallocate( JCoeff_L1 )
      deallocate( JCoeff_L2 );      deallocate( JCoeff_L3 )
! ---
    end subroutine set_porb_rot_matrix_jrep3

    subroutine set_porb_rot_matrix_jeff
      integer :: ia, il, im1, im2, im3, im4, im0, is1, is2, istmp
      integer :: immax_p, immax_d, lwork, info, size2
      complex(kind=CMPLDP) :: ztmp
      real(kind=DP), allocatable :: eigenval(:), rwork(:)
      complex(kind=CMPLDP), allocatable :: MatA(:,:), work(:), zwk(:,:)

      immax_p = 3;    immax_d = 5

      size2 = immax_p *ndim_spinor
      lwork = 2 *size2

      allocate( work(lwork) )
      allocate( eigenval(size2) ); eigenval = 0.0d0
      allocate( rwork( 3*size2 ) );
      allocate( MatA(size2,size2) )
      allocate( zwk(immax_p*ndim_spinor,immax_p*ndim_spinor) )

      Do ia=1, natm
         if (iproj_group(ia) == 0) cycle
         Do il=1, nloc
            if ( il /= 2+1 ) cycle
! --- t2g ----
            MatA = 0.0d0
            Do im1=1, immax_p
               Do im2=1, immax_p
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1-1 )*ndim_spinor +is2
                        ztmp = 0.0d0
                        Do im3=1, immax_d
                           Do im4=1, immax_d
                              ztmp = ztmp +porb_rot_matrix_real(ia,il,im3,im1+2) &
                             &      *Mat_LS_with_real_ylm_L2(im3,im4,istmp) &
                             &      *porb_rot_matrix_real(ia,il,im4,im2+2)
                           End Do
                        End Do
                        MatA( immax_p*(is1-1)+im1, immax_p*(is2-1)+im2 ) = ztmp
                     End do
                  End Do
               end Do
            End Do

            call zheev( 'V', 'U', size2, MatA, size2, &
                 &       eigenval, work, lwork, rwork, info )

            Do im0=1, immax_p*ndim_spinor
               Do im1=1, immax_d
                  Do is1=1, ndim_spinor
                     ztmp = 0.0d0
                     Do im2=1, immax_p
                        ztmp = ztmp +MatA(immax_p*(is1-1)+im2,im0) &
                             &      *porb_rot_matrix_real(ia,il,im1,im2+2)
                     End Do
                     porb_rot_matrix_cmplx(ia,il,immax_d*(is1-1)+im1,im0) &
                          & = ztmp
                  End Do
               End Do
            End Do
! --- eg ----
            Do im1=immax_p+1, immax_d
               Do is1=1, ndim_spinor
                  im0 = (im1-1)*ndim_spinor+is1
                  Do im2=1, immax_d
                     porb_rot_matrix_cmplx(ia,il,immax_d*(is1-1)+im2,im0 ) &
                          & = porb_rot_matrix_real(ia,il,im2,im1-immax_p)
                  End do
               End Do
            ENd Do

         End Do
      End Do
      deallocate( MatA ); deallocate( zwk )
      deallocate( work ); deallocate( rwork ); deallocate( eigenval )
! ---
    end subroutine set_porb_rot_matrix_jeff

    subroutine set_dm
      dm(:,:,:,:,1) =   dm_ssrep(:,:,:,:,1) +dm_ssrep(:,:,:,:,4)
      dm(:,:,:,:,2) =   dm_ssrep(:,:,:,:,2) +dm_ssrep(:,:,:,:,3)
      dm(:,:,:,:,3) =  (dm_ssrep(:,:,:,:,2) -dm_ssrep(:,:,:,:,3)) *zi
      dm(:,:,:,:,4) =   dm_ssrep(:,:,:,:,1) -dm_ssrep(:,:,:,:,4)
    end subroutine set_dm

    subroutine set_porb2_from_spin_densmat
      integer :: ia, it, il, im1, im2, im3, is1, is2, is3, istmp, jj, size2, immax
      complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )
      complex(kind=CMPLDP) :: ztmp, z1, z2
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:)

      call m_ES_set_Pauli_Matrix( PauliMatrix )

      do ia=1,natm
         it = ityp(ia)
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1
            size2 = immax *ndim_spinor
            allocate( zwk(size2,size2,ndim_magmom ) ); zwk = 0.0d0

            Do jj=1, ndim_magmom
               Do im1=1, immax
                  Do im2=1, immax
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           ztmp = 0.0d0
                           Do is3=1, ndim_spinor
                              istmp = ( is1 -1 )*ndim_spinor + is3
                              ztmp = ztmp +dm_ssrep(ia,il,im1,im2,istmp) &
                                   &      *PauliMatrix(jj,is3,is2)
                           End Do
                           zwk(im1+immax*(is1-1),im2+immax*(is2-1),jj) &
                                &  = ztmp
                        End Do
                     End Do
                  End Do
               ENd Do
            End Do
            Do jj=1, ndim_magmom
               Do im1=1, size2
                  ztmp  = 0.0d0
                  DO is2=1, ndim_spinor
                     DO is3=1, ndim_spinor
                        Do im2=1, immax
                           Do im3=1, immax
                              z1 = porb_rot_matrix_cmplx( ia,il,immax*(is2-1)+im2,im1 )
                              z2 = porb_rot_matrix_cmplx( ia,il,immax*(is3-1)+im3,im1 )
                              ztmp = ztmp + conjg(z1) &
                                   &       *zwk(immax*(is2-1)+im2,immax*(is3-1)+im3,jj) &
                                   &       *z2
                             End Do
                          End Do
                     End Do
                  End DO
                  porb2( ia, il, im1, jj ) = ztmp
               End Do
            End Do
            deallocate( zwk )
         End Do
      End do
    end subroutine set_porb2_from_spin_densmat

    subroutine set_porb2_from_chgmag_densmat
      integer :: ia, il, immax, size2, im1, im2, im3, is
      real(kind=DP) :: ctmp, cwk(ndim_chgpot)

      do ia=1,natm
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1
            size2 = immax

            Do im1=1, size2
               DO is=1, ndim_magmom
                  ctmp = 0.0d0
                  Do im2=1, immax
                     Do im3=1, immax
                        ctmp = ctmp &
                             &    + porb_rot_matrix_real( ia,il,im2,im1 )  &
                             &     *dm(ia,il,im2,im3,is) &
                             &     *porb_rot_matrix_real( ia,il,im3,im1 )
                     End Do
                  End Do
                  cwk(is) = ctmp
               End DO
               porb2( ia, il, im1, 1:ndim_magmom ) = cwk(1:ndim_magmom)
            End Do
         End Do
      ENd do

    end subroutine set_porb2_from_chgmag_densmat

    subroutine diagonalize_spin_densmat
      integer :: ia, il, im1, im2, is1, is2, lwork, size2, immax, istmp, info
      real(kind=DP), allocatable :: rwork(:), eigenval(:)
      complex(kind=CMPLDP), allocatable :: amat(:,:), work(:)
!
      do ia=1,natm
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1
            size2 = immax *ndim_spinor

            lwork = 2 *size2

            allocate( work(lwork) )
            allocate( amat(size2,size2) ); amat = 0.0d0
            allocate( eigenval(size2) ); eigenval = 0.0d0
            allocate( rwork( 3*size2 ) );

            amat = 0.0d0
            Do im1=1, immax
               Do im2=1, immax
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2
                        amat(im1+immax*(is1-1),im2+immax*(is2-1)) &
                             &  = dm_ssrep(ia,il,im1,im2,istmp)
                     End Do
                  End Do
               ENd Do
            End Do

            call zheev( 'V', 'U', size2, amat, size2, &
                 &       eigenval, work, lwork, rwork, info )

            if (info /= 0) then
               write(nfout,*) 'zheev : info=',info
            end if

            Do im1=1, size2
               Do im2=1, size2
!                  porb_rot_matrix_cmplx( ia,il,im2,im1 ) = amat(im2,im1) &
!                       &            * sqrt( eigenval(im1) )
                  porb_rot_matrix_cmplx( ia,il,im2,im1 ) = amat(im2,im1)
               ENd Do
            End Do
! --
            deallocate( amat );
            deallocate(rwork,work)
            deallocate( eigenval );
         end do
      end do
    end subroutine diagonalize_spin_densmat

    subroutine diagonalize_chg_densmat
      integer :: ia, il, immax, size2, lwork, im1, im2, info
      real(kind=DP) :: abstol,nfound, ctmp

      integer, allocatable :: iwork(:),ifail(:)
      real(kind=DP), allocatable :: rwork(:), amat(:,:), evec(:), rvmat(:,:)
      real(kind=DP), external :: dlamch
!
      do ia=1,natm
         if (iproj_group(ia) == 0) cycle

         Do il=1, lmax
            immax = 2*il -1;
            lwork = 8 *immax
            allocate( rwork(lwork) );
            allocate( iwork(5*immax) );   allocate( ifail(immax) )
            allocate( amat(immax,immax) )
            allocate( evec(immax) )
            allocate( rvmat(immax,immax) )

            amat = 0.0d0
            Do im1=1, immax
               Do im2=1, immax
                  amat(im1,im2) = dm(ia,il,im1,im2,1)
               ENd Do
            End Do
            call dsyevx( 'V', 'A', 'U', immax, amat, immax, 0.d0, 0.d0, 0, 0, &
                 &       abstol, nfound, evec, rvmat, immax, rwork, lwork, &
                 &       iwork, ifail, info )
            if (info /= 0) then
               write(nfout,*) 'dsyevx: info=',info
            end if

            Do im1=1, immax
               Do im2=1, immax
                  porb_rot_matrix_real( ia,il,im2,im1 ) = rvmat(im2,im1)
               ENd Do
            End Do
            deallocate(rwork);  deallocate(iwork);   deallocate(ifail)
            deallocate(amat);   deallocate(evec);    deallocate(rvmat)
         end do
      end do
    end subroutine diagonalize_chg_densmat

  end subroutine m_OP_rotation_noncl

  subroutine m_OP_alloc_compri_rot
    use m_Control_Parameters,  only : sw_diagonalize_population, population_diag_mode
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX
    integer :: num_orbitals, kfac

    kfac = 1
    if ( noncol ) then
       if ( sw_diagonalize_population == ON ) then
          if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX ) then
             kfac = ndim_spinor
          endif
       endif
    endif
    num_orbitals = nlmta_phi *kfac
    allocate(compr_rot_l(np_e,num_orbitals,nopr,ista_k:iend_k))
    allocate(compi_rot_l(np_e,num_orbitals,nopr,ista_k:iend_k))
    compr_rot_l = 0.d0; compi_rot_l = 0.0d0

  end subroutine m_OP_alloc_compri_rot

  subroutine m_OP_dealloc_compri_rot
    if (allocated(compr_rot_l) ) deallocate(compr_rot_l)
    if (allocated(compi_rot_l) ) deallocate(compi_rot_l)
  end subroutine m_OP_dealloc_compri_rot

  subroutine m_OP_set_compri_rot
    use m_Control_Parameters,  only : population_diag_mode
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX

    if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX ) then
       call set_compri_rot_chgden_mode
    else
       call set_compri_rot_spnden_mode
    endif
  end subroutine m_OP_set_compri_rot

  subroutine set_compri_rot_chgden_mode
    integer :: ik, ib, ia, it, immax
    integer :: lmt1, lmt2, il1, il2, im1, im2, tau1, tau2, iorb1, iorb2
    integer :: is1, is2, jj1, jj2, istmp
    complex(kind=CMPLDP), allocatable :: ztmp(:), zwork(:,:,:)

    allocate( ztmp(nopr) )
    allocate( zwork(nlmta_phi,nopr,1) );

    do ik = 1, kv3
       if(map_k(ik) /= myrank_k) cycle

       do ib = 1, np_e
          zwork = 0.0d0;
          do iorb1=1, nlmta_phi
             zwork(iorb1,:,1) = cmplx( compr_l(ib,iorb1,:,ik), &
                  &                    compi_l(ib,iorb1,:,ik) )
          end do
          do ia=1,natm
             it = ityp(ia)
             if (iproj_group(ia) == 0) cycle

             do lmt1=1,ilmt_phi(it)
                il1 = ltp_phi(lmt1,it);      im1 = mtp_phi(lmt1,it)
                tau1 = taup_phi(lmt1,it);    iorb1 = lmta_phi(lmt1,ia)
                immax = 2*il1 -1

                ztmp = 0.0d0
                do lmt2=1, ilmt_phi(it)
                   il2 = ltp_phi(lmt2,it);    im2 = mtp_phi(lmt2,it)
                   tau2 = taup_phi(lmt2,it);  iorb2 = lmta_phi(lmt2,ia)
                   if ( il1 /= il2 ) cycle
                   if ( tau1 /= tau2 ) cycle
                   ztmp(:) = ztmp(:) &
                        &    + porb_rot_matrix_real(ia,il1,im2,im1) &
                        &                 *zwork(iorb2,:,1)
                end do
                compr_rot_l(ib,iorb1,:,ik) = real(ztmp(:))
                compi_rot_l(ib,iorb1,:,ik) = aimag(ztmp(:))
             end do
          end do
       end do
    end do
    deallocate( ztmp );      deallocate( zwork )
  end subroutine set_compri_rot_chgden_mode

  subroutine set_compri_rot_spnden_mode
    integer :: ik, ib, ia, it, immax
    integer :: lmt1, lmt2, il1, il2, im1, im2, tau1, tau2, iorb1, iorb2
    integer :: is1, is2, jj1, jj2, istmp
    complex(kind=CMPLDP), allocatable :: ztmp(:,:), zwork(:,:,:)

    allocate( ztmp(nopr,ndim_spinor) )
    allocate( zwork(nlmta_phi,nopr,ndim_spinor) );

    do ik = 1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle

       do ib = 1, np_e
          zwork = 0.0d0;
          do iorb1=1, nlmta_phi
             zwork(iorb1,1:nopr,1) = cmplx( compr_l(ib,iorb1,1:nopr,ik), &
                  &                    compi_l(ib,iorb1,1:nopr,ik) )
             zwork(iorb1,1:nopr,2) = cmplx( compr_l(ib,iorb1,1:nopr,ik+1), &
                  &                    compi_l(ib,iorb1,1:nopr,ik+1) )
          end do

          do ia=1,natm
             it = ityp(ia)
             if (iproj_group(ia) == 0) cycle

             do lmt1=1,ilmt_phi(it)
                il1 = ltp_phi(lmt1,it);      im1 = mtp_phi(lmt1,it)
                tau1 = taup_phi(lmt1,it);    iorb1 = lmta_phi(lmt1,ia)
                immax = 2*il1 -1

                Do is1=1, ndim_spinor
                   ztmp = 0.0d0;

                   do lmt2=1, ilmt_phi(it)
                      il2 = ltp_phi(lmt2,it);    im2 = mtp_phi(lmt2,it)
                      tau2 = taup_phi(lmt2,it);  iorb2 = lmta_phi(lmt2,ia)
                      if ( il1 /= il2 ) cycle
                      if ( tau1 /= tau2 ) cycle

                      Do is2=1, ndim_spinor
                         jj2 = immax*(is2-1) +im2
                         jj1 = (im1 -1)*ndim_spinor +is1
                         ztmp(:,is2) = ztmp(:,is2) &
                              &      + conjg(porb_rot_matrix_cmplx(ia,il1,jj2,jj1)) &
                              &        *zwork(iorb2,:,is2)
                      end do
                   end do
                   compr_rot_l(ib,ndim_spinor*(iorb1-1)+is1,:,ik) &
                        &          = real(ztmp(:,1))
                   compi_rot_l(ib,ndim_spinor*(iorb1-1)+is1,:,ik) &
                        &          = aimag(ztmp(:,1))
                   compr_rot_l(ib,ndim_spinor*(iorb1-1)+is1,:,ik+1) &
                        &          = real(ztmp(:,2))
                   compi_rot_l(ib,ndim_spinor*(iorb1-1)+is1,:,ik+1) &
                       &          = aimag(ztmp(:,2))
                end do
             end do
          end do
       end do
    end do
    deallocate( ztmp );    deallocate( zwork )
  end subroutine set_compri_rot_spnden_mode

! ========================================================
  subroutine m_OP_init_Wfn_orb_proj(nfout)
    use m_Control_Parameters,  only : population_diag_mode, sw_diagonalize_population, &
         &                            sw_read_orb_rot_mat_file, &
         &                            sw_calc_score_sigma_bond, &
         &                            sw_write_rotated_orbitals
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX, LOCAL_POINT_GROUP
    use m_Files,   only : nfporb_dens_mat, nfporb_rot_mat, &
         &                m_Files_open_nfporb_dens_mat, m_Files_close_nfporb_dens_mat, &
         &                m_Files_open_nfporb_rot_mat, m_Files_close_nfporb_rot_mat

    integer, intent(in) :: nfout

    integer :: lmax, mmax

    lmax = nloc;   mmax = 2*lmax -1

    if ( sw_diagonalize_population == ON ) then
       if ( sw_read_orb_rot_mat_file == OFF ) then
          if ( noncol ) then
             call read_porb_dens_mat_ssrep
             call m_OP_rotation_noncl(nfout)
             deallocate(dm_ssrep)
          else
             call read_porb_dens_mat
             call m_OP_rotation_col(nfout)
             deallocate(dm)
          endif
       else
          if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX ) then
             call read_porb_rot_matrix_real
          else if ( population_diag_mode == LOCAL_POINT_GROUP ) then
             allocate( orb_irrep(natm,lmax,mmax) )
             call m_OP_calc_symm
          else
             if ( noncol ) call read_porb_rot_matrix_cmplx
          endif
          if ( sw_calc_score_sigma_bond == ON ) then
             if ( allocated(score_sigma_bond) ) deallocate(score_sigma_bond)
             allocate( score_sigma_bond(natm,lmax,mmax) );   score_sigma_bond = 0.0d0
             call m_OP_calc_score_sigma_bond
          endif
       endif
    endif
    if ( sw_write_rotated_orbitals == ON ) call m_OP_wd_phirt2_rotated

  contains

    subroutine read_porb_dens_mat
      integer :: ia, il, immax, im, im2, size2, is
      integer :: lmax, mmax, lun, num

      lmax = nloc;   mmax = 2*lmax -1

      if ( allocated(dm) ) deallocate(dm)
      allocate( dm(natm,lmax,mmax,mmax,nspin) ); dm = 0.0d0

      if ( mype == 0 ) then
         call m_Files_open_nfporb_dens_mat(1)
         lun = nfporb_dens_mat

         Do ia=1, natm
            if (iproj_group(ia) == 0) cycle
            read(lun) num
            Do is=1, nspin
               read(lun) num
               Do il=1, nloc
                  immax = 2*il -1
                  read(lun) num
                  Do im=1, immax
                     read(lun) ( dm(ia,il,im2,im,is), im2=1, immax )
                  End Do
               End do
            End Do
         End Do
         call m_Files_close_nfporb_dens_mat
      endif
      size2 = natm *nloc *mmax *mmax *nspin
      call mpi_bcast( dm, size2, mpi_double_precision, &
           &          0, MPI_CommGroup, ierr )

    end subroutine read_porb_dens_mat

    subroutine read_porb_dens_mat_ssrep
      integer :: ia, il, immax, im, im2, size2, is
      integer :: lmax, mmax, lun, num

      lmax = nloc;   mmax = 2*lmax -1

      if ( allocated(dm_ssrep) ) deallocate(dm_ssrep)
      allocate( dm_ssrep(natm,lmax,mmax,mmax,ndim_chgpot) ); dm_ssrep = 0.0d0

      if ( mype == 0 ) then
         call m_Files_open_nfporb_dens_mat(1)
         lun = nfporb_dens_mat

         Do ia=1, natm
            if (iproj_group(ia) == 0) cycle
            read(lun) num
            Do is=1, ndim_chgpot
               read(lun) num
               Do il=1, nloc
                  immax = 2*il -1
                  read(lun) num
                  Do im=1, immax
                     read(lun) ( dm_ssrep(ia,il,im2,im,is), im2=1, immax )
                  End Do
               End do
            End Do
         End Do
         call m_Files_close_nfporb_dens_mat
      endif
      size2 = natm *nloc *mmax *mmax *ndim_chgpot *2
      call mpi_bcast( dm_ssrep, size2, mpi_double_precision, &
           &          0, MPI_CommGroup, ierr )

    end subroutine read_porb_dens_mat_ssrep

    subroutine read_porb_rot_matrix_real
      integer :: ia, il, immax, im, im2, size2, lmax, mmax, num
      integer :: lun

      lmax = nloc;   mmax = 2*lmax -1

      if ( allocated(porb_rot_matrix_real) ) deallocate( porb_rot_matrix_real )
      allocate( porb_rot_matrix_real(natm,lmax,mmax,mmax) );
      porb_rot_matrix_real = 0.0d0

      if ( mype == 0 ) then
         call m_Files_open_nfporb_rot_mat(1)
         lun = nfporb_rot_mat
         read(lun) num

         Do ia=1, natm
            if (iproj_group(ia) == 0) cycle
            read(lun) num
            Do il=1, nloc
               immax = 2*il -1
               read(lun) num
               Do im=1, immax
                  read(lun) ( porb_rot_matrix_real(ia,il,im2,im), im2=1, immax )
               End do
            End Do
         End Do
         call m_Files_close_nfporb_rot_mat
      endif
      size2 = natm *nloc *mmax *mmax
      call mpi_bcast( porb_rot_matrix_real, size2, mpi_double_precision, &
           &          0, MPI_CommGroup, ierr )
!
    end subroutine read_porb_rot_matrix_real

    subroutine read_porb_rot_matrix_cmplx
      integer :: ia, il, immax, im, im2, size2, lmax, mmax, num, fac
      integer :: lun

      lmax = nloc;   mmax = 2*lmax -1
      fac = ndim_spinor

      if ( allocated(porb_rot_matrix_cmplx) ) deallocate( porb_rot_matrix_cmplx )
      allocate( porb_rot_matrix_cmplx(natm,lmax,fac*mmax,fac*mmax) );
      porb_rot_matrix_cmplx = 0.0d0

      if ( mype == 0 ) then
         call m_Files_open_nfporb_rot_mat(1)
         lun = nfporb_rot_mat
         read(lun) num

         Do ia=1, natm
            if (iproj_group(ia) == 0) cycle
            read(lun) num
            Do il=1, nloc
               immax = 2*il -1
               read(lun) num
               Do im=1, immax*fac
                  read(lun) ( porb_rot_matrix_cmplx(ia,il,im2,im), im2=1, immax*fac )
               End do
            End Do
         End Do
         call m_Files_close_nfporb_rot_mat
      endif

      size2 = natm *nloc *mmax *mmax *fac *fac *2
      call mpi_bcast( porb_rot_matrix_cmplx, size2, mpi_double_precision, &
           &          0, MPI_CommGroup, ierr )

    end subroutine read_porb_rot_matrix_cmplx

  end subroutine m_OP_init_Wfn_orb_proj

  subroutine m_OP_wd_Wfn_orb_proj
    use m_Control_Parameters,  only : population_diag_mode, sw_diagonalize_population, &
         &                            use_rotated_compri, SpinOrbit_mode, ekmode, &
         &                            wf_orb_proj_print_format, &
         &                            sw_calc_score_sigma_bond
    use m_Const_Parameters,    only : DIAG_CHARGE_DENSITY_MATRIX, &
         &                            DIAG_SPIN_DENSITY_MATRIX, &
         &                            DIAG_LS_with_t2g_octa, DIAG_LS, zi, Neglected, &
         &                            LOCAL_POINT_GROUP
    use m_PseudoPotential, only :    nloc, ntau, iproj_phi, lmtt_phi
    use m_Files,               only : m_Files_open_nfporb_rot_mat, &
         &                            m_Files_close_nfporb_rot_mat, nfporb_rot_mat, &
         &                            nfwfk_orb_proj, &
         &                            m_Files_open_nfwfk_orb_proj, &
         &                            m_Files_close_nfwfk_orb_proj

    integer :: lmax, mmax
    integer :: neg_t, kfac, num_orbitals
    real(kind=DP), allocatable :: compr(:,:,:,:), compi(:,:,:,:), norm_phig_mpi(:,:)

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX ) then
          if ( noncol ) kfac = ndim_spinor
       endif
    endif

    lmax = nloc;   mmax = ( 2*lmax -1 )*kfac
    num_orbitals = nlmta_phi *kfac

    if ( sw_diagonalize_population == ON .and. use_rotated_compri == YES ) then
       call m_OP_alloc_compri_rot
       call m_OP_set_compri_rot
    endif

    allocate(compr(neg,num_orbitals,1,kv3));  compr = 0.d0
    allocate(compi(neg,num_orbitals,1,kv3));  compi = 0.d0
    if(.not.allocated(norm_phig_mpi)) allocate(norm_phig_mpi(nlmtt_phi,kv3/nspin))
    norm_phig_mpi=0.d0

    call set_array_compri_etc

    neg_t = neg -num_extra_bands

    if ( mype == 0 ) then
       if ( ekmode == ON ) then
          if ( nk_in_the_process == 1 ) then
             call m_Files_open_nfwfk_orb_proj(2)
             call print_header
          else
             call m_Files_open_nfwfk_orb_proj(3)
          endif
       else
          call m_Files_open_nfwfk_orb_proj(2)
          call print_header
       endif
    endif

    if ( sw_diagonalize_population == ON ) then
       if ( use_rotated_compri == YES ) then
          if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX ) then
             call case_mode1
          else
             if ( noncol ) call case_mode2
          endif
       else
          if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX &
               & .or. population_diag_mode == LOCAL_POINT_GROUP ) then
             call case_mode1
          else
             if ( noncol ) call case_mode2
          endif
       endif
    else
       if ( SpinOrbit_Mode /= Neglected .and. wf_orb_proj_print_format == 1 ) then
          call case_with_j
       else
          call case_mode1
       endif
    endif

    if ( mype == 0 ) call m_Files_close_nfwfk_orb_proj

    deallocate( compr ); deallocate( compi ); deallocate( norm_phig_mpi )

    call m_OP_dealloc_compri_rot

  contains

    subroutine print_header
      integer :: num, iorb, ia, il, im, tau

      num = 0
      Do iorb=1, nlmta_phi
         call m_PP_tell_iorb_ia_l_m_tau(iorb,ia,il,im,tau)
         if ( iproj_group(ia) == 0 ) cycle
         num = num +1
      End Do

      write(nfwfk_orb_proj,'(A)') '# Orbital Projection for bands '
      write(nfwfk_orb_proj,*)
      write(nfwfk_orb_proj,'(A,I8)') 'num_kpoints = ', max(kv3,kv3_ek) /ndim_spinor
      write(nfwfk_orb_proj,'(A,I8)') 'num_bands   = ', neg_t
      write(nfwfk_orb_proj,'(A,I8)') 'nspin       = ', nspin /ndim_spinor
      write(nfwfk_orb_proj,'(A,I8)') 'num of orbitals = ', num *kfac
      write(nfwfk_orb_proj,*)
    end subroutine print_header

    subroutine set_array_compri_etc
      integer :: ik, ie, ib, iksnl
      integer :: iorb, lmt
      integer :: ia, il, im, tau, is
      real(kind=DP), allocatable :: compr_mpi(:,:,:,:), compi_mpi(:,:,:,:), &
           &                        norm_phig_mpi2(:,:)
      real(kind=DP), allocatable :: porb(:)

      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         iksnl = (ik-1)/nspin + 1
         if ( sw_diagonalize_population==ON .and. use_rotated_compri==YES ) then
!            do ie = ista_e, iend_e, istep_e
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               ib = map_z(ie)
               compr(ie,1:num_orbitals,1,ik) = compr_rot_l(ib,1:num_orbitals,1,ik)
               compi(ie,1:num_orbitals,1,ik) = compi_rot_l(ib,1:num_orbitals,1,ik)
             end do
         else
!            do ie = ista_e, iend_e, istep_e
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               ib = map_z(ie)
               compr(ie,1:num_orbitals,1,ik) = compr_l(ib,1:num_orbitals,1,ik)
               compi(ie,1:num_orbitals,1,ik) = compi_l(ib,1:num_orbitals,1,ik)
             end do
         endif
         norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
      end do

      if ( npes >1 ) then
         allocate( compr_mpi( neg, num_orbitals, 1, kv3 ) ); compr_mpi = 0.0d0
         allocate( compi_mpi( neg, num_orbitals, 1, kv3 ) ); compi_mpi = 0.0d0
         allocate( norm_phig_mpi2( nlmtt_phi, kv3/nspin ) )
         call mpi_allreduce( compr, compr_mpi, neg*num_orbitals*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( compi, compi_mpi, neg*num_orbitals*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*kv3/nspin, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         compr = compr_mpi /dble(nrank_g);   compi = compi_mpi/dble(nrank_g)
         norm_phig_mpi = norm_phig_mpi2 /dble(nrank_e) /dble(nrank_g)
         deallocate( compr_mpi ); deallocate( compi_mpi );
         deallocate( norm_phig_mpi2 )
      end if
    end subroutine set_array_compri_etc

! ==========
    subroutine case_mode1
      integer :: ik, iorb, ia, il, im, tau, ib
      integer :: iksnl, lmtt, it, is

      real(kind=DP), allocatable :: porb(:,:)
      logical, save :: First = .true.

      if ( mype /= 0 ) return

      if ( First ) then
         write(nfwfk_orb_proj,'(A,I3)') "population_diag_mode = ", population_diag_mode
         write(nfwfk_orb_proj,'(A)')
         write(nfwfk_orb_proj,'(A)') "# Orbital Info."
         if ( sw_diagonalize_population == ON ) then
            write(nfwfk_orb_proj,'(A)',advance='no') &
                 &    "  iorb  ia  l  m' tau  element  key"
            if ( sw_calc_score_sigma_bond == ON ) then
               write(nfwfk_orb_proj,'(A)',advance='no') "  score_bond"
            endif
            if ( population_diag_mode == LOCAL_POINT_GROUP ) then
               write(nfwfk_orb_proj,'(A)') "  symmetry"
            else
               write(nfwfk_orb_proj,*)
            endif
         else
            write(nfwfk_orb_proj,'(A)') "  iorb  ia  l  m  tau  element  key"
         endif

         Do iorb=1, num_orbitals
            call m_PP_tell_iorb_ia_l_m_tau(iorb,ia,il,im,tau)
            if ( iproj_group(ia) == 0) cycle
            it = ityp(ia)
            if ( sw_diagonalize_population == ON ) then
               write(nfwfk_orb_proj,'(I5,I5,I3,I3,2X,I3,4X,A5,I4)',advance='no') &
                    &       iorb, ia, il-1, im, tau, speciesname(it), atom_key(ia)
               if ( sw_calc_score_sigma_bond == ON ) then
                  write(nfwfk_orb_proj,'(1X,F10.5)',advance='no') &
                    &       score_sigma_bond(ia,il,im)
               endif
               if ( population_diag_mode == LOCAL_POINT_GROUP ) then
                  write(nfwfk_orb_proj,'(5x,a8)') orb_irrep(ia,il,im)
               else
                  write(nfwfk_orb_proj,*)
               endif
            else
               write(nfwfk_orb_proj,'(I5,I5,I3,I3,I3,4X,A5,I6)') &
                    &       iorb, ia, il-1, im, tau, speciesname(it), atom_key(ia)
            endif
         End Do
         write(nfwfk_orb_proj,*)
         First = .false.
      endif

      allocate( porb( neg,nlmta_phi) )

      do ik = 1, kv3, ndim_spinor
         if ( ekmode == ON ) then
            if ( ik+nk_in_the_process -1 > kv3_ek ) cycle
         endif

         iksnl = (ik-1)/nspin + 1

         write(nfwfk_orb_proj,'(A)') "================= "

         if ( sw_band_unfolding == ON ) then
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz_refcell(ik+nk_in_the_process-1,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik, " ( ", vkxyz_refcell(ik,1:3,BUCS), " )"
            endif
         else
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz(ik,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik, " ( ", vkxyz(ik,1:3,BUCS), " )"
            endif
         endif

         porb = 0.0d0
         if ( sw_diagonalize_population == OFF .or. use_rotated_compri == YES ) then
            do iorb = 1,nlmta_phi
               call m_PP_tell_iorb_lmtt(iorb,lmtt)
               Do is=1, ndim_spinor
                  do ib = 1, neg
                     porb(ib,iorb) = porb(ib,iorb) &
                       & + ( compr(ib,iorb,1,ik+is-1)**2 &
                       &    +compi(ib,iorb,1,ik+is-1)**2 ) &
                       &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                  end do
               end do
            End Do
         else
            call set_porb_mode1( ik, porb )
         endif

         do iorb = 1,nlmta_phi
            call m_PP_tell_iorb_ia_l_m_tau(iorb,ia,il,im,tau)
            if ( iproj_group(ia) == 0) cycle

            if ( sw_diagonalize_population == ON ) then
               write(nfwfk_orb_proj,'(I5,I5,3I3,A)') iorb, ia, il-1, im, tau, &
                    &                             " : iorb, ia, l, m', tau"
            else
               write(nfwfk_orb_proj,'(I5,I5,3I3,A)') iorb, ia, il-1, im, tau, &
                    &                             " : iorb, ia, l, m, tau"
            endif
            write(nfwfk_orb_proj,'(4F18.10)') ( porb(ib,iorb ), ib=1, neg_t )
         end do
         write(nfwfk_orb_proj,*)
      end do
      deallocate( porb )

    end subroutine case_mode1

    subroutine set_porb_mode1( ik, porb )
      integer, intent(in) :: ik
      real(kind=DP), intent(out) :: porb( neg,nlmta_phi*ndim_spinor )

      integer :: ia, ib, it, lmax, mmax, immax, ismax
      integer :: lmt, il, im, tau, ip, iorb, lmtt
      integer :: lmt1, lmt2, im1, im2, iorb1, iorb2, lmtt1, lmtt2
      integer :: is, is1, is2, istmp, iksnl, im3
      real(kind=DP) :: c1, c2, ctmp
      complex(kind=CMPLDP) :: ztmp, z1, z2

      complex(kind=CMPLDP), allocatable :: dm_ssrep(:,:,:,:,:)
      real(kind=DP), allocatable :: porb0(:,:,:,:), dm(:,:,:,:,:)

      lmax = nloc;   mmax = 2*lmax -1
      iksnl = (ik-1)/nspin + 1

      if ( noncol ) then
         allocate( dm_ssrep( ntau, lmax, mmax, mmax, ndim_chgpot ) );
         allocate( dm( ntau, lmax, mmax, mmax, ndim_magmom ) );
         allocate( porb0( ntau, lmax, mmax, ndim_magmom ) );
      else
         allocate( dm_ssrep( ntau, lmax, mmax, mmax, 1 ) );
         allocate( dm( ntau, lmax, mmax, mmax, 1 ) );
         allocate( porb0( ntau, lmax, mmax, 1 ) );
      endif

      ibloop: do ib=1, neg
         do ia=1,natm
            it = ityp(ia)
            if(iproj_group(ia) == 0) cycle

            dm_ssrep = 0.0d0;  dm = 0.0d0;  porb0 = 0.d00

            ! diagonal part
            do lmt=1,ilmt_phi(it)
               il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
               tau = taup_phi(lmt,it)
               ip  = iproj_phi(lmt,it)
               iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

               if ( noncol ) then
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2

                        z1 = dcmplx( compr(ib,iorb,1,ik+is1-1 ), &
                             &       compi(ib,iorb,1,ik+is1-1 ) )
                        z2 = dcmplx( compr(ib,iorb,1,ik+is2-1 ), &
                             &       compi(ib,iorb,1,ik+is2-1 ) )
                        ztmp = z1 *conjg(z2) &
                             &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                        dm_ssrep(tau,il,im,im,istmp) = dm_ssrep(tau,il,im,im,istmp) &
                             &                        +ztmp
                     end do
                  end Do
               else
                  z1 = dcmplx( compr(ib,iorb,1,ik ), &
                       &       compi(ib,iorb,1,ik ) )
                  z2 = dcmplx( compr(ib,iorb,1,ik ), &
                       &       compi(ib,iorb,1,ik ) )
                  ztmp = z1 *conjg(z2) &
                       &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                  dm_ssrep(tau,il,im,im,1) = dm_ssrep(tau,il,im,im,1) +ztmp
               endif
            end do
            ! non-diagonal part
            do lmt2=1,ilmt_phi(it)
               do lmt1=1,ilmt_phi(it)
                  if ( lmt1 == lmt2 ) cycle
                  if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
                  if (taup_phi(lmt1,it) /= taup_phi(lmt2,it)) cycle

                  il = ltp_phi(lmt1,it);     tau = taup_phi(lmt1,it)
                  ip  = iproj_phi(lmt1,it)
                  im1   = mtp_phi(lmt1,it);     im2 = mtp_phi(lmt2,it);
                  iorb1 = lmta_phi(lmt1,ia);  iorb2 = lmta_phi(lmt2,ia)
                  lmtt1 = lmtt_phi(lmt1,it);  lmtt2 = lmtt_phi(lmt2,it)

                  c1 = ( qorb(iorb1) +qorb(iorb2) ) /2.0d0      ! approx
                  c2 = sqrt(norm_phig_mpi(lmtt1,iksnl)*norm_phig_mpi(lmtt2,iksnl))

                  if ( noncol ) then
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           istmp = ( is1 -1 )*ndim_spinor + is2
                           z1 = dcmplx( compr(ib,iorb1,1,ik+is1-1 ), &
                                &       compi(ib,iorb1,1,ik+is1-1 ) )
                           z2 = dcmplx( compr(ib,iorb2,1,ik+is2-1 ), &
                                &       compi(ib,iorb2,1,ik+is2-1 ) )
                           ztmp = z1 *conjg(z2) *( 1.0d0 +c1 /c2 )
                           dm_ssrep(tau,il,im1,im2,istmp) &
                                &  = dm_ssrep(tau,il,im1,im2,istmp) + ztmp
                        end do
                     end do
                  else
                     z1 = dcmplx( compr(ib,iorb1,1,ik ), &
                          &       compi(ib,iorb1,1,ik ) )
                     z2 = dcmplx( compr(ib,iorb2,1,ik ), &
                          &       compi(ib,iorb2,1,ik ) )
                     ztmp = z1 *conjg(z2) *( 1.0d0 +c1 /c2 )
                     dm_ssrep(tau,il,im1,im2,1) = dm_ssrep(tau,il,im1,im2,1) +ztmp
                  endif
               end do
            end do

            if ( noncol ) then
               dm(:,:,:,:,1) =  dm_ssrep(:,:,:,:,1) +dm_ssrep(:,:,:,:,4)
               dm(:,:,:,:,2) =  dm_ssrep(:,:,:,:,2) +dm_ssrep(:,:,:,:,3)
               dm(:,:,:,:,3) = (dm_ssrep(:,:,:,:,2) -dm_ssrep(:,:,:,:,3)) *zi
               dm(:,:,:,:,4) =  dm_ssrep(:,:,:,:,1) -dm_ssrep(:,:,:,:,4)
               ismax = 1
            else
               dm(:,:,:,:,1) = dm_ssrep(:,:,:,:,1)
               ismax = 1
            endif
            ! ----
            Do tau=1, ntau
               Do il=1, lmax
                  immax = 2 *il -1
                  Do im1=1, immax
                     DO is=1, ismax
                        ctmp = 0.0d0
                        Do im2=1, immax
                           Do im3=1, immax
                              ctmp = ctmp &
                                   &    +porb_rot_matrix_real( ia,il,im2,im1 )  &
                                   &       *dm(tau,il,im2,im3,is) &
                                   &       *porb_rot_matrix_real( ia,il,im3,im1 )
                           End Do
                        End Do
                        porb0( tau, il, im1, is ) = ctmp
                     End DO
                  End DO
               End Do
            End Do

            Do lmt=1,ilmt_phi(it)
               il = ltp_phi(lmt,it);  im  = mtp_phi(lmt,it);  iorb = lmta_phi(lmt,ia)
               tau = taup_phi(lmt,it)
               porb( ib,iorb ) = porb0( tau, il,im,1 )
            ENd do
         End do
      End do ibloop
      deallocate( dm ); deallocate( dm_ssrep ); deallocate( porb0 )

    end subroutine set_porb_mode1
! ==========

    subroutine case_mode2
      integer :: ik, iorb1, ia, il, im, tau, iksnl, lmtt
      integer :: im0, iorb0, my_l, ib, it
      real(kind=DP) :: val_j, val_mj
      complex(kind=CMPLDP) :: z1, z2

      real(kind=DP), allocatable :: porb(:,:)
      logical, save :: First = .true.

      if ( mype /= 0 ) return

      if ( First ) then
         write(nfwfk_orb_proj,'(A,I3)') "population_diag_mode = ", population_diag_mode
         write(nfwfk_orb_proj,'(A)')
         write(nfwfk_orb_proj,'(A)') "# Orbital Info."

         if ( population_diag_mode == DIAG_LS ) then
            write(nfwfk_orb_proj,'(A)') &
                 &                 "  iorb  ia  l  ms'  tau elment   j    mj"

            do iorb1 = 1, nlmta_phi *ndim_spinor
               iorb0 = int( (iorb1-1)/ndim_spinor )+1
               call m_PP_tell_iorb_ia_l_m_tau(iorb0,ia,il,im0,tau)
               if ( iproj_group(ia) == 0) cycle
               im = ( im0 -1 )*ndim_spinor +mod( iorb1 -1, ndim_spinor ) +1
               it = ityp(ia);  my_l = il -1

               if ( my_l > 0 ) then
                  if ( im <= 2*my_l ) then
                     val_j = my_l -0.5d0;    val_mj = -val_j +(im -1)
                  else
                     val_j = my_l +0.5d0;    val_mj = -val_j +(im -2*my_l -1)
                  endif
               else
                  val_j = my_l +0.5d0;       val_mj = -val_j +(im -2*my_l -1)
               endif
               write(nfwfk_orb_proj,'(I5,I5,I3,I3,I3,4X,A5,X,F4.1,X,F4.1)') &
               &       iorb1, ia, il-1, im, tau, speciesname(it), &
               &       val_j, val_mj
            end do

         else
            write(nfwfk_orb_proj,'(A)',advance='no') &
                 &           "  iorb  ia  l  ms'   tau  element  key"
            if ( sw_calc_score_sigma_bond == ON ) then
               write(nfwfk_orb_proj,'(A)',advance='no') "  score_bond"
            endif
            if ( population_diag_mode == LOCAL_POINT_GROUP ) then
               write(nfwfk_orb_proj,'(A)') "  symmetry"
            else
               write(nfwfk_orb_proj,*)
            endif

            do iorb1 = 1, nlmta_phi *ndim_spinor
               iorb0 = int( (iorb1-1)/ndim_spinor )+1
               call m_PP_tell_iorb_ia_l_m_tau(iorb0,ia,il,im0,tau)
               if ( iproj_group(ia) == 0) cycle
               im = ( im0 -1 )*ndim_spinor +mod( iorb1 -1, ndim_spinor ) +1
               it = ityp(ia)
               write(nfwfk_orb_proj,'(I5,I5,I3,I3,I3,4X,A5,I6)') &
                    &       iorb1, ia, il-1, im, tau, speciesname(it), atom_key(ia)
               if ( sw_calc_score_sigma_bond == ON ) then
                  write(nfwfk_orb_proj,'(3X,F10.5)',advance='no') &
                    &       score_sigma_bond(ia,il,im)
               endif
               if ( population_diag_mode == LOCAL_POINT_GROUP ) then
                  write(nfwfk_orb_proj,'(5x,a8)') orb_irrep(ia,il,im)
               else
                  write(nfwfk_orb_proj,*)
               endif
            end do
         endif
         write(nfwfk_orb_proj,*)
         First = .false.
      endif

      allocate( porb( neg,nlmta_phi*ndim_spinor ) )

      do ik = 1, kv3, ndim_spinor
         if ( ekmode == ON ) then
            if ( ik+nk_in_the_process -1 > kv3_ek ) cycle
         endif

         iksnl = (ik-1)/nspin + 1

         write(nfwfk_orb_proj,'(A)') "================= "
         if ( sw_band_unfolding == ON ) then
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz_refcell(ik+nk_in_the_process-1,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik, " ( ", vkxyz_refcell(ik,1:3,BUCS), " )"
            endif
         else
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz(ik,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik, " ( ", vkxyz(ik,1:3,BUCS), " )"
            endif
         endif

         porb = 0.0d0
         if ( use_rotated_compri == YES ) then
            do iorb1 = 1, nlmta_phi *ndim_spinor
               iorb0 = int( (iorb1-1)/ndim_spinor )+1
               call m_PP_tell_iorb_lmtt(iorb0,lmtt)
               do ib = 1, neg
                  z1 = dcmplx( compr(ib,iorb1,1,ik ),   compi(ib,iorb1,1,ik   ) )
                  z2 = dcmplx( compr(ib,iorb1,1,ik+1),  compi(ib,iorb1,1,ik+1 ) )
                  porb(ib,iorb1) = porb(ib,iorb1) &
                       &   + (z1+z2)*conjg(z1+z2) &
                       &     *( 1.d0+qorb(iorb0)/norm_phig_mpi(lmtt,iksnl) )
               end do
            end do
         else
            call set_porb_mode2( ik, porb )
         endif

         do iorb1 = 1, nlmta_phi *ndim_spinor
            iorb0 = int( (iorb1-1)/kfac )+1
            call m_PP_tell_iorb_ia_l_m_tau(iorb0,ia,il,im0,tau)
            if ( iproj_group(ia) == 0) cycle

            im = ( im0 -1 )*ndim_spinor +mod( iorb1 -1, ndim_spinor ) +1

            if ( population_diag_mode == DIAG_LS ) then
               my_l = il -1
               if ( my_l > 0 ) then
                  if ( im <= 2*my_l ) then
                     val_j = my_l -0.5d0;     val_mj = -val_j +(im -1)
                  else
                     val_j = my_l +0.5d0;     val_mj = -val_j +(im -2*my_l -1)
                  endif
               else
                  val_j = my_l +0.5d0;        val_mj = -val_j +(im -2*my_l -1)
               endif
               write(nfwfk_orb_proj,'(I5,I5, 3I3,A,A,F4.1,A,F4.1,A)') &
                    &                iorb1, ia, il-1, im, tau, &
                    &                ' : iorb, ia, l, ms'', tau', &
                    &                "     ( j= ", val_j, "  mj= ", val_mj, " )"
            else
               write(nfwfk_orb_proj,'(I5,I5,3I3,A)') &
                    &                iorb1, ia, il-1, im, tau, &
                    &                ' : iorb, ia, l, ms'', tau'
            endif
            write(nfwfk_orb_proj,'(4F18.10)') &
                 &              ( porb(ib,iorb1 ), ib=1, neg_t )
         End Do
         write(nfwfk_orb_proj,*)
      end do
      deallocate( porb )

    end subroutine case_mode2

    subroutine set_porb_mode2( ik, porb )
      integer, intent(in) :: ik
      real(kind=DP), intent(out) :: porb( neg,nlmta_phi*ndim_spinor )

      integer :: ia, ib, it, lmax, mmax, immax, ismax
      integer :: lmt, il, im, tau, ip, iorb, lmtt
      integer :: lmt1, lmt2, im1, im2, iorb1, iorb2, lmtt1, lmtt2
      integer :: is, is1, is2, is3, istmp, iksnl, im3
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: ztmp, z1, z2, z3
      real(kind=DP), allocatable :: porb0(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: dm_ssrep(:,:,:,:,:), dm(:,:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:)
      complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

      lmax = nloc;   mmax = 2*lmax -1

      allocate( dm_ssrep( ntau, lmax, mmax, mmax, ndim_chgpot ) );
      allocate( porb0( ntau, lmax, mmax*ndim_spinor, ndim_magmom ) );
      allocate( zwk( mmax*ndim_spinor, mmax*ndim_spinor ) );

      call m_ES_set_Pauli_Matrix( PauliMatrix )

      iksnl = (ik-1)/nspin + 1

      ibloop: do ib=1, neg
         do ia=1,natm
            it = ityp(ia)
            if(iproj_group(ia) == 0) cycle

            dm_ssrep = 0.0d0;  dm = 0.0d0;  porb0 = 0.d00

            ! diagonal part
            do lmt=1,ilmt_phi(it)
               il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
               tau = taup_phi(lmt,it)
               ip  = iproj_phi(lmt,it)
               iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     istmp = ( is1 -1 )*ndim_spinor + is2

                     z1 = dcmplx( compr(ib,iorb,1,ik+is1-1 ), &
                          &       compi(ib,iorb,1,ik+is1-1 ) )
                     z2 = dcmplx( compr(ib,iorb,1,ik+is2-1 ), &
                          &       compi(ib,iorb,1,ik+is2-1 ) )
                     ztmp = z1 *conjg(z2) &
                          &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                     dm_ssrep(tau,il,im,im,istmp) = dm_ssrep(tau,il,im,im,istmp) &
                          &                        +ztmp
                  end do
               end Do
            end do
            ! non-diagonal part
            do lmt2=1,ilmt_phi(it)
               do lmt1=1,ilmt_phi(it)
                  if ( lmt1 == lmt2 ) cycle
                  if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
                  if (taup_phi(lmt1,it) /= taup_phi(lmt2,it)) cycle

                  il = ltp_phi(lmt1,it);     tau = taup_phi(lmt1,it)
                  ip  = iproj_phi(lmt1,it)
                  im1   = mtp_phi(lmt1,it);     im2 = mtp_phi(lmt2,it);
                  iorb1 = lmta_phi(lmt1,ia);  iorb2 = lmta_phi(lmt2,ia)
                  lmtt1 = lmtt_phi(lmt1,it);  lmtt2 = lmtt_phi(lmt2,it)

                  c1 = ( qorb(iorb1) +qorb(iorb2) ) /2.0d0      ! approx
                  c2 = sqrt(norm_phig_mpi(lmtt1,iksnl)*norm_phig_mpi(lmtt2,iksnl))

                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2
                        z1 = dcmplx( compr(ib,iorb1,1,ik+is1-1 ), &
                             &       compi(ib,iorb1,1,ik+is1-1 ) )
                        z2 = dcmplx( compr(ib,iorb2,1,ik+is2-1 ), &
                             &       compi(ib,iorb2,1,ik+is2-1 ) )
                        ztmp = z1 *conjg(z2) *( 1.0d0 +c1 /c2 )
                        dm_ssrep(tau,il,im1,im2,istmp) &
                             &  = dm_ssrep(tau,il,im1,im2,istmp) + ztmp
                     end do
                  end do
               end do
            end do

            ismax = 1
            ! ----
            Do tau=1, ntau
               Do il=1, lmax
                  immax = 2 *il -1
                  Do is=1, ndim_magmom
                     zwk = 0.0d0

                     Do im1=1, immax
                        Do im2=1, immax
                           Do is1=1, ndim_spinor
                              Do is2=1, ndim_spinor
                                 ztmp = 0.0d0
                                 Do is3=1, ndim_spinor
                                    istmp = ( is1 -1 )*ndim_spinor + is3
                                    ztmp = ztmp +dm_ssrep(tau,il,im1,im2,istmp) &
                                         &      *PauliMatrix(is,is3,is2)
                                 End Do
                                 zwk(im1+immax*(is1-1),im2+immax*(is2-1)) = ztmp
                              End Do
                           End Do
                        End Do
                     ENd Do
                     Do im1=1, immax *ndim_spinor
                        ztmp = 0.0d0
                        DO is2=1, ndim_spinor
                           DO is3=1, ndim_spinor
                              Do im2=1, immax
                                 Do im3=1, immax
                                    z1 = porb_rot_matrix_cmplx(ia,il,immax*(is2-1)+im2,im1)
                                    z2 = zwk(immax*(is2-1)+im2,immax*(is3-1)+im3)
                                    z3 = porb_rot_matrix_cmplx(ia,il,immax*(is3-1)+im3,im1)
                                    ztmp = ztmp +conjg(z1) *z2 *z3
                                 End Do
                              End Do
                           End Do
                        End DO
                        porb0( tau, il, im1, is ) = ztmp
                     End Do
                  End Do
               End Do
            End Do

            Do lmt=1,ilmt_phi(it)
               il = ltp_phi(lmt,it);  im1 = mtp_phi(lmt,it);  iorb1 = lmta_phi(lmt,ia)
               tau = taup_phi(lmt,it)
               Do is1=1, ndim_spinor
                  im2 = (im1 -1)*ndim_spinor +is1
                  iorb2 = (iorb1 -1)*ndim_spinor +is1
                  porb( ib,iorb2 ) = porb0( tau, il,im2,1 )
               ENd do
            End Do
         End do
      End do ibloop
      deallocate( dm_ssrep ); deallocate( zwk ); deallocate( porb0 )

    end subroutine set_porb_mode2
! =========

    subroutine case_with_j
      integer :: ik, iksnl, ia, ig, ii, it, ilp, ll, tau, ip
      integer :: iorb, lmtt1, ib, m1
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2

      real(kind=DP), allocatable :: porb(:)
      complex(kind=CMPLDP), allocatable :: zcomp(:,:,:)

      logical, save :: First = .true.

      if ( mype /= 0 ) return

      allocate( porb(neg) )

      Do ik=1, kv3, ndim_spinor
         if ( ekmode == ON ) then
            if ( ik+nk_in_the_process -1 > kv3_ek ) cycle
         endif

         write(nfwfk_orb_proj,'(A)') "================= "

         if ( sw_band_unfolding == ON ) then
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz_refcell(ik+nk_in_the_process-1,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                 &            ik, " ( ", vkxyz_refcell(ik,1:3,BUCS), " )"
            endif
         else
            if ( ekmode == ON ) then
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &            ik +nk_in_the_process -1, " ( ", &
                    &            vkxyz(ik,1:3,BUCS), " )"
            else
               write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
                 &            ik, " ( ", vkxyz(ik,1:3,BUCS), " )"
            endif
         endif

         iksnl = ( ik -1 )/nspin +1

         Do ia=1, natm
            ig = iproj_group(ia)
            if ( ig == 0 ) cycle

            do ii=1,num_proj_elems(ig)
               ip = proj_group( ii, ig )
               it = proj_attribute(ip)%ityp
               ilp = proj_attribute(ip)%l +1
               ll = proj_attribute(ip)%l
               tau = proj_attribute(ip)%t
!
               allocate( zcomp( -ll:ll, neg, ndim_spinor ) ); zcomp = 0.0d0
               call tranform_compri_r2c_sph( ia, it, ll, tau, ik, &
                    &                        compr, compi, zcomp )

               if ( ll == 0 ) then
                  call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                  call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                  write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                       &         ia, ll+0.5d0, ll, 0.5d0, tau, ' : ia, j, l, mj, tau'
                  Do ib=1, neg
                     z1 = cmplx( compr(ib,iorb,1,ik),   compi(ib,iorb,1,ik) )
                     z2 = cmplx( compr(ib,iorb,1,ik+1), compi(ib,iorb,1,ik+1) )
                     porb(ib) = ( z1*conjg(z1) +z2*conjg(z2) ) &
                          &     *( 1.d0+qorb(iorb)/norm_phig_mpi( lmtt1, iksnl ) )
                  End Do
                  write(nfwfk_orb_proj,'(4F18.10)') &
                       &               ( porb( neordr(ib,ik) ), ib=1, neg_t )
               else
! j_up
                  Do m1=-ll -1, ll
                     c1 = dble( ll +m1 + 1 ) / dble( 2 *ll +1 )
                     c2 = dble( ll -m1 )     / dble( 2 *ll +1 )
                     c1 = sqrt(c1);  c2 = sqrt(c2)

                     call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                     call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                     write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                          &         ia, ll+0.5d0, ll, m1+0.5d0, tau, &
                          &             ' : ia, j, l, mj, tau'
                     Do ib=1, neg
                        z1 = 0.0d0
                        if ( m1 > -ll -1 ) z1 = z1 +c1 *zcomp( m1,    ib, 1 )
                        if ( m1 < ll     ) z1 = z1 +c2 *zcomp( m1 +1, ib, 2 )
                        porb(ib) = z1 *conjg(z1) &
                             &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt1,iksnl) )
                     End Do
                     write(nfwfk_orb_proj,'(4F18.10)') &
                          &         ( porb( neordr(ib,ik) ), ib=1, neg_t )
                  End Do
! j_down
                  Do m1=-ll+1, ll
                     c1 = dble( ll -m1 + 1 ) / dble( 2 *ll +1 )
                     c2 = dble( ll +m1 )     / dble( 2 *ll +1 )
                     c1 = sqrt(c1);  c2 = -sqrt(c2)

                     call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                     call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                     write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                          &         ia, ll-0.5d0, ll, m1-0.5d0, tau, &
                          &             ' : ia, j, l, mj, tau'
                     Do ib=1, neg
                        z1 = c1 *zcomp( m1-1, ib, 1 ) +c2 *zcomp( m1, ib, 2 )
                        porb(ib) = z1 *conjg(z1) &
                             &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt1,iksnl) )
                     End Do
                     write(nfwfk_orb_proj,'(4F18.10)') &
                          &         ( porb(ib), ib=1, neg_t )
                  End Do
               end if
               deallocate( zcomp )

            End Do
         End Do
         write(nfwfk_orb_proj,*)

      End Do
      deallocate( porb )

    end subroutine case_with_j

    subroutine tranform_compri_r2c_sph( ia, it, ll, tau, ik, &
         &                              compr, compi, zcomp )
      integer, intent(in) :: ia, it, ll, tau, ik
      real(kind=DP), intent(in) :: compr( neg, nlmta_phi, 1, kv3 )
      real(kind=DP), intent(in) :: compi( neg, nlmta_phi, 1, kv3 )
      complex(kind=CMPLDP), intent(out) :: zcomp( -ll:ll, neg, ndim_spinor )

      integer :: m1, m2, ib, iorb
      complex(kind=CMPLDP) :: z1
      complex(kind=CMPLDP) :: ztmp( ndim_spinor )

      Do m1=-ll, ll
         Do m2=1, 2*ll +1
            call find_iorb_from_lmt( ia, it, ll, m2, tau, iorb )
            Do ib=1, neg
               ztmp(1) = cmplx( compr(ib,iorb,1,ik),   compi(ib,iorb,1,ik)   )
               ztmp(2) = cmplx( compr(ib,iorb,1,ik+1), compi(ib,iorb,1,ik+1) )

               if ( ll == 0 ) z1 = MatU_ylm_RC_L0( m2, m1 )
               if ( ll == 1 ) z1 = MatU_ylm_RC_L1( m2, m1 )
               if ( ll == 2 ) z1 = MatU_ylm_RC_L2( m2, m1 )
               if ( ll == 3 ) z1 = MatU_ylm_RC_L3( m2, m1 )

               zcomp(m1,ib,:) = zcomp(m1,ib,:) +z1 *ztmp(:)
            End Do
         End Do
      End Do
    end subroutine tranform_compri_r2c_sph

    subroutine find_iorb_from_lmt( ia, it, ll, mm, tau, iorb )
      integer, intent(in) :: ia, it, ll, mm, tau
      integer, intent(out) :: iorb

      integer :: lmt1, l1, m1, t1

      iorb = 0
      Do lmt1=1, ilmt_phi(it)
         l1 = ltp_phi(lmt1,it); m1 = mtp_phi(lmt1,it);  t1 = taup_phi(lmt1,it)
         if ( l1 == ll +1 .and. m1 == mm .and. t1 == tau ) then
            exit
         endif
      ENd Do
      iorb = lmta_phi( lmt1,ia )

    end subroutine find_iorb_from_lmt

  end subroutine m_OP_wd_Wfn_orb_proj

  subroutine m_OP_wd_phirt2_rotated
    use m_PseudoPotential,  only : nlmt_phi, ilmt_phi, ltp_phi, mtp_phi, taup_phi, &
         &                         phirt, xh, radr_paw, rmax, nmesh
    use m_Ionic_System,  only : cps, pos
    use m_Control_Parameters,  only : population_diag_mode
    use m_Const_Parameters,   only : DIAG_CHARGE_DENSITY_MATRIX, &
         &                           DIAG_SPIN_DENSITY_MATRIX, DIAG_LS_with_t2g_octa, &
         &                           LOCAL_POINT_GROUP, LOCAL_DOUBLE_POINT_GROUP
!
    call change_of_coordinate_system(altv,pos,natm,natm,cps)

    select case( population_diag_mode )
    case ( DIAG_CHARGE_DENSITY_MATRIX )
       call case_spinor1_basis
    case ( DIAG_SPIN_DENSITY_MATRIX )
       call case_spinor2_basis
    case ( DIAG_LS_with_t2g_octa )
       call case_spinor2_basis
    case ( LOCAL_POINT_GROUP )
       call case_spinor1_basis
    case ( LOCAL_DOUBLE_POINT_GROUP )
       call case_spinor2_basis
    end select

  contains

    subroutine case_spinor1_basis
      integer :: ia, it
      integer :: lmt1, il1, im1, tau1, lmt2, il2, im2, tau2
      real(kind=DP) :: rcut, dx, dy, dz
      real(kind=DP) :: origin(3), vec_span(3,3)
      real(kind=DP), allocatable :: phi_on_mesh(:,:,:,:)
      real(kind=DP), allocatable :: phi2_on_mesh(:,:,:)

      real(kind=DP), allocatable :: phi_on_mesh_rot(:,:,:)

      integer :: lun = 4000
      integer :: nxmax, nymax, nzmax

      nxmax = 10;     nymax = 10;     nzmax = 10

      Do ia=1, natm
         it = ityp(ia)
         if (iproj_group(ia) == 0) cycle

         rcut = 3.0d0
         dx = rcut /nxmax;        dy = rcut /nymax;        dz = rcut /nzmax

         origin(1:3) = cps(ia,1:3) -rcut
         vec_span = 0.0d0
         vec_span(1,1) = rcut *2;   vec_span(2,2) = rcut *2;   vec_span(3,3) = rcut *2

         allocate( phi_on_mesh(nlmt_phi,-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax) )
         allocate( phi_on_mesh_rot(-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax) )
         allocate( phi2_on_mesh(-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax) )

         call set_phirt_on_mesh( it, -nxmax, nxmax, -nymax, nymax, -nzmax, nzmax, &
              &                  dx, dy, dz, phi_on_mesh )
         ! -
         do lmt1=1,ilmt_phi(it)
            il1 = ltp_phi(lmt1,it);  im1 = mtp_phi(lmt1,it);  tau1 = taup_phi(lmt1,it)

            phi_on_mesh_rot = 0.0d0

            Do lmt2=1, ilmt_phi(it)
               il2 = ltp_phi(lmt2,it);  im2 = mtp_phi(lmt2,it);  tau2 = taup_phi(lmt2,it)
               if ( il1 /= il2 ) cycle
               if ( tau1 /= tau2 ) cycle
               phi_on_mesh_rot(:,:,:) = phi_on_mesh_rot(:,:,:) &
                    &                  +porb_rot_matrix_real(ia,il1,im2,im1) &
                    &                   *phi_on_mesh(lmt2,:,:,:)
            End Do
            phi2_on_mesh(:,:,:) = phi_on_mesh_rot(:,:,:)  &
                 &               *phi_on_mesh_rot(:,:,:)

            if ( mype == 0 ) then
               call printout_phi2( lun, ia, il1, im1, tau1, &
                    &             -nxmax, nxmax, -nymax, nymax, -nzmax, nzmax, &
                    &              origin, vec_span, phi2_on_mesh, 1 )
            endif
         End do
         deallocate( phi_on_mesh );       deallocate( phi_on_mesh_rot )
         deallocate( phi2_on_mesh )
      End Do

    end subroutine case_spinor1_basis

    subroutine case_spinor2_basis
      integer :: ia, it
      integer :: lmt1, il1, im1, tau1, lmt2, il2, im2, tau2
      integer :: j1, j2, is1, is2, immax
      real(kind=DP) :: rcut, dx, dy, dz
      real(kind=DP) :: origin(3), vec_span(3,3)
      real(kind=DP), allocatable :: phi_on_mesh(:,:,:,:)
      real(kind=DP), allocatable :: phi2_on_mesh(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: phi_on_mesh_rot(:,:,:,:)

      integer :: lun = 4000
      integer :: nxmax, nymax, nzmax

      nxmax = 10;     nymax = 10;     nzmax = 10

      Do ia=1, natm
         it = ityp(ia)
         if (iproj_group(ia) == 0) cycle

         rcut = 3.0d0
         dx = rcut /nxmax;        dy = rcut /nymax;        dz = rcut /nzmax

         origin(1:3) = cps(ia,1:3) -rcut
         vec_span = 0.0d0
         vec_span(1,1) = rcut *2;   vec_span(2,2) = rcut *2;   vec_span(3,3) = rcut *2

         allocate( phi_on_mesh(nlmt_phi,-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax) )
         allocate( phi_on_mesh_rot(-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax,ndim_spinor) )
         allocate( phi2_on_mesh(-nxmax:nxmax,-nymax:nymax,-nzmax:nzmax,ndim_spinor) )

         call set_phirt_on_mesh( it, -nxmax, nxmax, -nymax, nymax, -nzmax, nzmax, &
              &                  dx, dy, dz, phi_on_mesh )
         ! -
         do lmt1=1,ilmt_phi(it)
            il1 = ltp_phi(lmt1,it);  im1 = mtp_phi(lmt1,it);  tau1 = taup_phi(lmt1,it)

            if ( population_diag_mode == DIAG_LS_with_t2g_octa ) then
               if ( il1 /= 2 +1 ) cycle
            endif

            Do is1=1, ndim_spinor
               phi_on_mesh_rot = 0.0d0

               Do lmt2=1, ilmt_phi(it)
                  il2 = ltp_phi(lmt2,it);  im2 = mtp_phi(lmt2,it);
                  tau2 = taup_phi(lmt2,it)
                  immax = 2 *il2 -1
                  if ( il1 /= il2 ) cycle
                  if ( tau1 /= tau2 ) cycle
                  Do is2=1, ndim_spinor
                     j1 = ( im1 -1 )*ndim_spinor +is1
                     j2 = immax *(is2 -1) +im2

                     phi_on_mesh_rot(:,:,:,is2) = phi_on_mesh_rot(:,:,:,is2) &
                       &                        +porb_rot_matrix_cmplx(ia,il1,j2,j1)&
                       &                        *phi_on_mesh(lmt2,:,:,:)
                  End do
               End Do
#if 0
               phi2_on_mesh(:,:,:,1) = conjg( phi_on_mesh_rot(:,:,:,1) ) &
                    &                  *phi_on_mesh_rot(:,:,:,1)
               phi2_on_mesh(:,:,:,2) = conjg( phi_on_mesh_rot(:,:,:,2) ) &
                    &                  *phi_on_mesh_rot(:,:,:,2)
               if ( mype == 0 ) then
                  call printout_phi2( lun, ia, il1, j1, tau1, &
                       &             -nxmax, nxmax, -nymax, nymax, -nzmax, nzmax, &
                       &              origin, vec_span, phi2_on_mesh, 2 )
               endif
#else
               phi2_on_mesh(:,:,:,1) = conjg( phi_on_mesh_rot(:,:,:,1) ) &
                    &                  *phi_on_mesh_rot(:,:,:,1) &
                    &                +conjg( phi_on_mesh_rot(:,:,:,2) ) &
!                    &                -conjg( phi_on_mesh_rot(:,:,:,2) ) &
                    &                  *phi_on_mesh_rot(:,:,:,2)
               if ( mype == 0 ) then
                  call printout_phi2( lun, ia, il1, j1, tau1, &
                       &             -nxmax, nxmax, -nymax, nymax, -nzmax, nzmax, &
                       &              origin, vec_span, phi2_on_mesh, 1 )
               endif
#endif
            End Do
         End do
         deallocate( phi_on_mesh );       deallocate( phi_on_mesh_rot )
         deallocate( phi2_on_mesh )
      End Do

    end subroutine case_spinor2_basis

    subroutine printout_phi2( lun, ia, il, im, tau, &
         &                    nxs, nxe, nys, nye, nzs, nze, &
         &                    origin, vec_span, phi2_on_mesh, ncomp )
      integer, intent(in) :: lun, ia, il, im, tau, ncomp
      integer, intent(in) :: nxs, nxe, nys, nye, nzs, nze
      real(kind=DP), intent(in) :: phi2_on_mesh(nxs:nxe,nys:nye,nzs:nze,ncomp)
      real(kind=DP), intent(in) :: origin(3), vec_span(3,3)

      integer :: i, nx, ny, nz
      character*4 char1
      character*1 char2, char4
      character*2, char3
      character*64 file1

      write(char1,'(I4.4)') ia
      write(char2,'(I1.1)') il -1
      write(char3,'(I2.2)') im
      write(char4,'(I1.1)') tau

      file1 = "phirot_squared.ia_" // char1 // ".l_" // char2 // ".tau_" &
           &                 // char4 // ".no_" // char3 // ".xsf"
      open( unit=lun, file=file1, status="unknown", form="formatted" )

      write(lun,'(A)') 'CRYSTAL'
      write(lun,'(A)') 'PRIMVEC'
      write(lun,'(3F20.10)') altv(1:3,1) *Bohr
      write(lun,'(3F20.10)') altv(1:3,2) *Bohr
      write(lun,'(3F20.10)') altv(1:3,3) *Bohr
      write(lun,'(A)') 'CONVVEC'
      write(lun,'(3F20.10)') altv(1:3,1) *Bohr
      write(lun,'(3F20.10)') altv(1:3,2) *Bohr
      write(lun,'(3F20.10)') altv(1:3,3) *Bohr
      write(lun,'(A)') 'PRIMCOORD'
      write(lun,*) natm, 1
      Do i=1, natm
         write(lun,'(I5,3F20.10)') nint(iatomn(ityp(i))), cps(i,1:3)*Bohr
      End Do
      write(lun,'(A)') "BEGIN_BLOCK_DATAGRID_3D"
      write(lun,'(A,I6,A,I4,A,I4,A,I4,A)') &
           &    "  phi2 on mesh ( ia=", ia, ", il=", il, ", tau=", tau, &
           &                      ", no=", im,  " )"
      write(lun,'(A)') "BEGIN_DATAGRID_3D_1"
      write(lun,*) nxe -nxs +1, nye -nys +1, nze -nzs +1
      write(lun,'(3F20.10)') origin*Bohr
      write(lun,'(3F20.10)') vec_span(1:3,1) *Bohr
      write(lun,'(3F20.10)') vec_span(1:3,2) *Bohr
      write(lun,'(3F20.10)') vec_span(1:3,3) *Bohr
      write(lun,'(5F10.5)') (((phi2_on_mesh(nx,ny,nz,1),nx=nxs,nxe), &
           &                  ny=nys,nye),nz=nzs,nze)
      write(lun,'(A)') "END_DATAGRID_3D"
      if ( ncomp == 2 ) then
        write(lun,'(A)') "BEGIN_DATAGRID_3D_2"
         write(lun,*) nxe -nxs +1, nye -nys +1, nze -nzs +1
         write(lun,'(3F20.10)') origin*Bohr
         write(lun,'(3F20.10)') vec_span(1:3,1) *Bohr
         write(lun,'(3F20.10)') vec_span(1:3,2) *Bohr
         write(lun,'(3F20.10)') vec_span(1:3,3) *Bohr
         write(lun,'(5F10.5)') (((phi2_on_mesh(nx,ny,nz,2),nx=nxs,nxe), &
              &                  ny=nys,nye),nz=nzs,nze)
         write(lun,'(A)') "END_DATAGRID_3D"
      endif
      write(lun,'(A)') "END_BLOCK_DATAGRID_3D"
      close(lun)
    end subroutine printout_phi2

    subroutine set_phirt_on_mesh( it, nxs, nxe, nys, nye, nzs, nze, &
         &                        dx, dy, dz, phi_on_mesh )
      integer, intent(in) :: it, nxs, nxe, nys, nye, nzs, nze
      real(kind=DP), intent(in) :: dx, dy, dz
      real(kind=DP), intent(out) :: phi_on_mesh(nlmt_phi,nxs:nxe,nys:nye,nzs:nze)

      integer :: lmt1, il1, im1, tau1, nx, ny, nz, ind, nspher1
      real(kind=DP) :: cx, cy, cz, dist, ylm1, c1, ctmp, f1, f2
      real(kind=DP) :: dist_min = 1.0D-10

      do lmt1=1,ilmt_phi(it)
         il1 = ltp_phi(lmt1,it); im1 = mtp_phi(lmt1,it); tau1 = taup_phi(lmt1,it)

         Do nz=nzs, nze
            Do ny=nys, nye
               Do nx=nxs, nxe
                  cx = nx *dx;     cy = ny *dy;       cz = nz *dz
                  dist = sqrt(cx**2 +cy**2 +cz**2)

                  if ( dist < dist_min ) then
                     dist = dist_min
                     ind = 1
                     c1 = phirt(ind,il1,tau1,it) /radr_paw(ind,it)
                  else
                     ctmp = log( dist /rmax(it) ) *xh(it) +nmesh(it)
                     ind = int(ctmp)
                     if ( ind < 1 ) then
                        ind = 1
                        c1 = phirt(ind,il1,tau1,it) /radr_paw(ind,it)
                     else
                        f1 = ctmp -int(ctmp)
                        f2 = 1.0d0 -f1
                        c1 = ( phirt(ind,il1,tau1,it)/radr_paw(ind,it) )**f2 &
                             & *( phirt(ind+1,il1,tau1,it)/radr_paw(ind+1,it)) **f1
                     endif
                  endif
                  ! -- ylm --
                  nspher1 = ( il1 -1 )**2 +im1
                  call sphr( 1, nspher1, cx, cy, cz, ylm1 )
                  ! -------
                  phi_on_mesh(lmt1,nx,ny,nz) = c1 *ylm1
               End do
            ENd Do
         End Do
      End do
    end subroutine set_phirt_on_mesh

  end subroutine m_OP_wd_phirt2_rotated

  subroutine m_OP_calc_score_sigma_bond
    use m_Control_Parameters, only : sw_diagonalize_population, population_diag_mode
    use m_Const_Parameters,   only : DIAG_CHARGE_DENSITY_MATRIX, &
         &                           LOCAL_POINT_GROUP

    integer :: ia, ja, nx, ny, nz
    integer :: lmax, mmax, immax
    integer :: il, l1, im1, im2, ig, ilm1, im3, is2, kfac
    integer :: nsph_min, natom_nn
    real(kind=DP) :: vec(3), dist, dist_min, factor, c1, rcut
    complex(kind=CMPLDP) :: z1

    real(kind=DP), allocatable :: val_sph(:), work(:,:)

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
            &  .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
          if ( noncol ) kfac = ndim_spinor
       endif
    endif

    lmax = nloc;   mmax = 2*lmax -1
    allocate( val_sph(mmax) )
    allocate( work(lmax,mmax*kfac) )

    factor = 1.2D0

    Do ia=1, natm
       ig = iproj_group(ia)
       if ( ig == 0 ) cycle

       dist_min = 1.0D2
       Loop_ja1: Do ja=1, natm
          Do nx=-1, 1
             Do ny=-1, 1
                Do nz=-1, 1
                   vec(1) = cps(ja,1) -cps(ia,1) &
                        &   +nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                   vec(2) = cps(ja,2) -cps(ia,2) &
                        &   +nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                   vec(3) = cps(ja,3) -cps(ia,3) &
                        &   +nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                   dist = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )
                   if ( dist > 0.01 ) dist_min = min( dist, dist_min )
                ENd Do
             End Do
          End Do
       End Do Loop_ja1
       rcut = dist_min *factor
!       write(*,*) "rcut = ", ia, rcut

       natom_nn = 0;        work = 0.0d0

       Loop_ja2: Do ja=1, natm
          Do nx=-1, 1
             Do ny=-1, 1
                Do nz=-1, 1
                   vec(1) = cps(ja,1) -cps(ia,1) &
                        &   +nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                   vec(2) = cps(ja,2) -cps(ia,2) &
                        &   +nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                   vec(3) = cps(ja,3) -cps(ia,3) &
                        &   +nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                   dist = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )

                   if ( dist > 0.01 .and. dist < rcut ) then
                      vec = vec /dist
                      natom_nn = natom_nn +1

                      Do il=1, lmax
                         l1 = il -1
                         nsph_min = l1**2 +1;
                         immax = 2*l1 +1

                         Do im1=1, immax
                            ilm1 = nsph_min +im1 -1
                            call sphr( 1, ilm1, vec(1), vec(2), vec(3), val_sph(im1) )
                         End Do

                         if ( sw_diagonalize_population == ON ) then
                            if ( population_diag_mode &
                                 &     == DIAG_CHARGE_DENSITY_MATRIX &
                                 & .or. population_diag_mode &
                                 &        == LOCAL_POINT_GROUP ) then
                               Do im1=1, immax
                                  c1 = 0.0
                                  Do im2=1, immax
                                     c1 = c1 +val_sph(im2) &
                                          &  *porb_rot_matrix_real(ia,il,im2,im1)
                                  End do
                                  work(il,im1) = work(il,im1) +c1**2 /immax
                               End Do
                            else
                               Do im1=1, immax *kfac
                                  z1 = 0.0
                                  Do is2=1, ndim_spinor
                                     Do im2=1, immax
                                        im3=immax*(is2-1) +im2
                                        z1 = z1 +val_sph(im2) &
                                             &  *porb_rot_matrix_cmplx(ia,il,im3,im1)
                                     End Do
                                  End do
                                  work(il,im1) = work(il,im1) +z1*conjg(z1) /immax
                               End Do
                            endif
                         else
                            Do im1=1, immax
                               work(il,im1) = work(il,im1) +val_sph(im1)**2 /immax
                            ENd Do
                         endif
                      End Do
                   end if
                End Do
             ENd Do
          End Do

       End Do Loop_ja2
       score_sigma_bond(ia,:,:) = work(:,:) *PAI4
! -----------------------
    End Do
#if 0
!    stop
    write(*,*) "UUUU"
    Do ia=1, natm
       Do il=1, lmax
          Do im1=1, mmax
             write(*,*) ia, il, im1, score_sigma_bond(ia,il,im1)
          ENd Do
       ENd Do
    End Do
#endif
    deallocate(val_sph); deallocate(work)

  end subroutine m_OP_calc_score_sigma_bond

  subroutine m_OP_calc_symm
    use m_Control_Parameters, only : sw_diagonalize_population, population_diag_mode, &
         &                           ndim_spinor
    use m_Const_Parameters,   only : DIAG_CHARGE_DENSITY_MATRIX, PAI, zi, &
         &                           LOCAL_DOUBLE_POINT_GROUP
    use m_Crystal_Structure,   only : il
    use m_CS_SpaceGroup,        only : m_CS_gen_opr_rspace_full, get_point_group_name
    use m_CS_Magnetic,    only : m_CS_gen_opr_spinor_full
    use m_Files,  only : nfout
    use m_Representation,  only :  get_irr_rep, get_irr_rep_dblegrp
    use m_SpinOrbit_Potential,  only : m_SO_set_MatU_ylm_RC, &
         &                             m_SO_calc_MatLS_orb_s_to_f,&
         &                             m_SO_diagonalize_MatLS
    use m_Ionic_System,  only : ionic_mass

    integer, parameter :: max_nn = 8
    character(len=9) :: system
    character(len=3) :: pg_name

    integer :: nsym, nsym0, natom_nn
    integer :: lmax, mmax, kfac
    real(kind=DP) :: factor, rcut, dist_min

    integer :: ia, ig, ja

    real(kind=DP) :: sym_axis(3,3), angle(3)
    integer, allocatable :: Jtype(:)
    real(kind=DP), allocatable :: Jmass(:), BondVec(:,:)

    real(kind=DP), allocatable, target :: rot_space_cubic(:,:,:), rot_space_hex(:,:,:)
    real(kind=DP), pointer :: rot_space_org(:,:,:)
    real(kind=DP), allocatable :: rot_space_wk(:,:,:)
    complex(kind=CMPLDP), allocatable, target :: &
         &                     rot_spinor_cubic(:,:,:), rot_spinor_hex(:,:,:)
    complex(kind=CMPLDP), allocatable :: rot_spinor_wk(:,:,:)
    complex(kind=CMPLDP), pointer :: rot_spinor_org(:,:,:)

    logical :: flg_dblegrp

    flg_dblegrp = .false.
    if ( noncol .and. population_diag_mode == LOCAL_DOUBLE_POINT_GROUP ) then
       flg_dblegrp = .true.
    endif

    allocate( rot_space_cubic(3,3,48) );   allocate( rot_space_hex(3,3,24) );
    if ( flg_dblegrp ) then
       allocate( rot_spinor_cubic(2,2,48) );  allocate( rot_spinor_hex(2,2,48) )
    endif

    call m_SO_set_MatU_ylm_RC
    if ( flg_dblegrp ) call set_Jcoeff_L0_L3

    system = 'cubic';  nsym = 48
    call m_CS_gen_opr_rspace_full( system, nsym, rot_space_cubic )
    if ( flg_dblegrp ) call m_CS_gen_opr_spinor_full( system, nsym, rot_spinor_cubic )

    system = 'hexagonal';  nsym = 24
    call m_CS_gen_opr_rspace_full( system, nsym, rot_space_hex )
    if ( flg_dblegrp ) call m_CS_gen_opr_spinor_full( system, nsym, rot_spinor_hex )

    allocate( BondVec(max_nn,3) );   allocate( Jtype(max_nn) )
    allocate( Jmass(max_nn) )

    lmax = nloc;   mmax = 2*lmax -1;     nsym0 = nsym

    factor = 1.2D0

    Do ia=1, natm
       ig = iproj_group(ia)
       if ( ig == 0 ) cycle

       nsym = nsym0
       call calc_dist_min( ia, dist_min )
       rcut = dist_min *factor

       call set_atoms_nn( ia, rcut, max_nn, natom_nn, BondVec, Jtype )
       call check_if_hexagonal( max_nn, natom_nn, BondVec, system )

       call set_symmetry_axis( max_nn, natom_nn, BondVec, sym_axis )
       call calc_euler_angles2( sym_axis, angle(1), angle(2), angle(3) )

       if ( system == 'cubic' ) then
          nsym0 = 48;  rot_space_org => rot_space_cubic
          if ( flg_dblegrp ) rot_spinor_org => rot_spinor_cubic
       else
          nsym0 = 24;  rot_space_org => rot_space_hex
          if ( flg_dblegrp ) rot_spinor_org => rot_spinor_hex
       endif

       allocate( rot_space_wk(3,3,nsym0) )
       if ( flg_dblegrp ) allocate( rot_spinor_wk(2,2,nsym0) )

       if ( ipriorb_rot > 2 ) then
          write(nfout,'(A,I5)') "** Symmetry Axis of Atom :", ia
          write(nfout,'(A,3F20.10)')  'vec1 : ',sym_axis(1:3,1)
          write(nfout,'(A,3F20.10)')  'vec2 : ',sym_axis(1:3,2)
          write(nfout,'(A,3F20.10)')  'vec3 : ',sym_axis(1:3,3)
          write(nfout,'(A)') "** Euler Angles (deg.)"
          write(nfout,'(A,3F20.8)') 'alpha beta gamma :', angle(1:3) /PAI *180.0d0
          write(nfout,'(A,A)') "** local system is ", system
       endif

       if ( flg_dblegrp ) then
          call set_porb_rot_matrix_cmplx0( ia )
          call calc_local_point_group( max_nn, natom_nn, sym_axis, BondVec, Jtype, &
               &                       nsym0, nsym, pg_name, &
               &                       rot_space_org, rot_spinor_org )
          call get_orb_irrep_dblegrp( ia, pg_name )
       else
          call set_porb_rot_matrix_real0( ia )
          call calc_local_point_group( max_nn, natom_nn, sym_axis, BondVec, Jtype, &
               &                       nsym0, nsym, pg_name, &
               &                       rot_space_org )
          call get_orb_irrep( ia, pg_name )
       endif
       deallocate( rot_space_wk )
       if ( flg_dblegrp ) deallocate( rot_spinor_wk )
    End Do
    !
    deallocate( BondVec )
    deallocate( Jtype );   deallocate( Jmass )
    deallocate( rot_space_cubic );   deallocate( rot_space_hex )
    if ( flg_dblegrp ) then
       deallocate( rot_spinor_cubic );   deallocate( rot_spinor_hex )
    endif

  contains

    subroutine check_if_hexagonal( max_nn, natom_nn, BondVec, system )
      integer, intent(in) :: max_nn, natom_nn
      real(kind=DP), intent(in) :: BondVec(max_nn,3)
      character(len=9), intent(out) :: system

      integer :: i, j, k
      real(kind=DP) :: vec(3), norm

      system = 'cubic'
      Do i=1, natom_nn-2
         Do j=i+1, natom_nn-1
            Do k=j+1, natom_nn
               vec(:) = BondVec(i,:) +BondVec(j,:) +BondVec(k,:)
               norm = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )
               if ( norm < 1.0E-3 ) system = 'hexagonal'
            End Do
         End Do
      ENd Do
    end subroutine check_if_hexagonal

    subroutine get_orb_irrep( ia, pg_name )
      integer, intent(in) :: ia
      character(len=3), intent(in) :: pg_name

      real(kind=DP), allocatable :: crotylm(:,:,:) ! d(mmax,nsph,nopr)
      integer, allocatable :: iylm(:,:,:) ! d(mmax,nsph,nopr)
      integer, allocatable :: nylm(:,:) ! d(nsph,nopr)
      real(kind=DP), allocatable :: matc(:,:,:), matwk(:,:)
      real(kind=DP), allocatable :: basis_vec(:,:), porb_rot_tmp(:,:)

      integer :: rep(48,12), nrep, nsym, nd
      integer :: num_basis, nadd
      character(len=3) :: symbol_irrep(12)
      character(len=4) :: active_rep(12)

      integer :: il1, l1, im1, im2, im3, im4, immax
      integer :: iopr, isph1, isph2, nsph_min, irep, iy, irep_cand
      real(kind=DP) :: c1, fac, ctmp, norm

      if ( pg_name == "Dih" .or. pg_name == "Civ" ) then
         orb_irrep(ia,:,:) = pg_name // ' ' // 'NON '
         return
      else
         call get_irr_rep( pg_name, nsym, nrep, rep, symbol_irrep, active_rep )
      endif

      if(.not.allocated(crotylm)) allocate(crotylm(mmax,lmax**2,nsym))
      if(.not.allocated(iylm))    allocate(iylm(mmax,lmax**2,nsym))
      if(.not.allocated(nylm))    allocate(nylm(lmax**2,nsym))

      crotylm = 0.0d0
      call get_crotylm( lmax, mmax, lmax**2, nsym, crotylm, iylm, nylm, rot_space_wk )

      Do il1=1, lmax
         l1 = il1 -1
         nsph_min = l1**2 +1;
         immax = 2*l1 +1

         allocate( matc(immax,immax,nsym) );  matc = 0.0d0

         Do iopr=1, nsym
            Do im1=1, immax
               isph1 = nsph_min +im1 -1
               Do iy=1, nylm(isph1,iopr)
                  isph2 = iylm(iy,isph1,iopr)
                  im2 = isph2 -nsph_min +1
                  matc(im2,im1,iopr) = crotylm(iy,isph1,iopr)
               End Do
            End Do
         ENd Do

         num_basis = 0
         allocate( basis_vec(immax,immax) ); basis_vec = 0.0d0
         allocate( matwk(immax,immax) )

         Do irep=1, nrep
            nd = rep(1,irep)
            matwk = 0.0d0
            Do im1=1, immax
               Do im2=1, immax
                  ctmp = 0.0d0
                  Do iopr=1, nsym
                     fac = rep(iopr,irep)
                     c1 = matc(im2,im1,iopr)
                     ctmp = ctmp +fac *c1
                  End Do
                  norm = ctmp *dble(nd) /dble(nsym)
                  if ( abs(norm) > 1.0D-6 ) matwk(im2,im1) = norm
               End Do
            End Do
            call count_basis( immax, matwk, num_basis, basis_vec, nadd )
            Do im1=1, nadd
               orb_irrep(ia,il1,num_basis+im1) = pg_name // ' ' // symbol_irrep(irep)
            End Do
            num_basis = num_basis +nadd
         End Do
         if ( num_basis /= immax ) call phase_error_with_msg(nfout,"num_basis /= immax",__LINE__,__FILE__)

         allocate( porb_rot_tmp(immax,immax) )
         porb_rot_tmp(1:immax,1:immax) &
              &     = porb_rot_matrix_real(ia,il1,1:immax,1:immax)

         Do im1=1, immax
            Do im2=1, immax
               ctmp = 0.0d0
               Do im3=1, immax
                  ctmp = ctmp +basis_vec(im3,im1) *porb_rot_tmp(im2,im3)
               End Do
               porb_rot_matrix_real(ia,il1,im2,im1) = ctmp
            End Do
         End Do

         deallocate( matwk );      deallocate( matc )
         deallocate( basis_vec );  deallocate( porb_rot_tmp )
      End Do
      deallocate(crotylm);     deallocate(iylm);     deallocate(nylm)

    end subroutine get_orb_irrep

    subroutine get_orb_irrep_dblegrp( ia, pg_name )
      integer, intent(in) :: ia
      character(len=3), intent(in) :: pg_name

      real(kind=DP), allocatable :: crotylm(:,:,:) ! d(mmax,nsph,nopr)
      integer, allocatable :: iylm(:,:,:) ! d(mmax,nsph,nopr)
      integer, allocatable :: nylm(:,:) ! d(nsph,nopr)
      real(kind=DP), allocatable :: basis_vec(:,:), matwk(:,:)
      complex(kind=CMPLDP), allocatable :: matc(:,:,:), porb_rot_tmp(:,:)
      complex(kind=CMPLDP), allocatable :: zsum(:), work(:,:)

      integer :: nrep, nsym, nd, iloop, istart, iend
      integer :: num_basis, nadd

      complex(kind=CMPLDP) :: rep_dblegrp(48,2,16)
      character(len=4) :: symbol_irrep_dblegrp(16)

      integer :: il1, l1, im1, im2, im3, im4, immax, is1, is2
      integer :: iopr, isph1, isph2, nsph_min, irep, iy, irep_cand
      real(kind=DP) :: norm, max_norm
      complex(kind=CMPLDP) :: z1, ztmp, fac

      if ( pg_name == "Dih" .or. pg_name == "Civ" ) then
         orb_irrep(ia,:,:) = pg_name // ' ' // 'NON '
         return
      else
         call get_irr_rep_dblegrp( pg_name, nsym, nrep, rep_dblegrp, &
              &                    symbol_irrep_dblegrp )
      endif

      if(.not.allocated(crotylm)) allocate(crotylm(mmax,lmax**2,nsym))
      if(.not.allocated(iylm))    allocate(iylm(mmax,lmax**2,nsym))
      if(.not.allocated(nylm))    allocate(nylm(lmax**2,nsym))

      crotylm = 0.0d0
      call get_crotylm( lmax, mmax, lmax**2, nsym, crotylm, iylm, nylm, rot_space_wk )

      Do il1=1, lmax
         l1 = il1 -1
         nsph_min = l1**2 +1;
         immax = 2*l1 +1

         allocate( zsum(immax*ndim_spinor) )
         allocate( matc(immax*ndim_spinor,immax*ndim_spinor,nsym) );  matc = 0.0d0
         allocate( work(immax*ndim_spinor,immax*ndim_spinor) )

         Do iopr=1, nsym
            work = 0.0d0
            Do im1=1, immax
               isph1 = nsph_min +im1 -1
               Do iy=1, nylm(isph1,iopr)
                  isph2 = iylm(iy,isph1,iopr)
                  im2 = isph2 -nsph_min +1
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        work( immax*(is2-1)+im2,immax*(is1-1)+im1 ) &
                             &    = crotylm(iy,isph1,iopr) *rot_spinor_wk(is2,is1,iopr)
                     End Do
                  ENd Do
               End Do
            End Do
            Do im1=1, immax*ndim_spinor
               Do im2=1, immax*ndim_spinor
                  ztmp = 0.0d0
                  Do im3=1, immax*ndim_spinor
                     Do im4=1, immax*ndim_spinor
                        select case (l1)
                        case (0)
                           ztmp = ztmp +work(im3,im4) *JCoeff_L0(im4,im2) &
                                &      *conjg(JCoeff_L0(im3,im1) )
                        case (1)
                           ztmp = ztmp +work(im3,im4) *JCoeff_L1(im4,im2) &
                           &           *conjg(JCoeff_L1(im3,im1) )
                        case (2)
                           ztmp = ztmp +work(im3,im4) *JCoeff_L2(im4,im2) &
                           &           *conjg(JCoeff_L2(im3,im1) )
                        case (3)
                           ztmp = ztmp +work(im3,im4) *JCoeff_L3(im4,im2) &
                           &           *conjg(JCoeff_L3(im3,im1) )
                        end select
                     End Do
                  End Do
                  matc( im1, im2, iopr ) = ztmp
               ENd Do
            End Do
         ENd Do

         num_basis = 0
         allocate( basis_vec(immax*ndim_spinor,immax*ndim_spinor) ); basis_vec = 0.0d0
         allocate( matwk(immax*ndim_spinor,immax*ndim_spinor) )

         Do iloop=1, 2
            if ( iloop == 1 ) then
               istart = 1;   iend = l1 *ndim_spinor
            else
               istart = l1 *ndim_spinor +1;  iend = immax *ndim_spinor
            endif

            Do irep=1, nrep
               nd = rep_dblegrp(1,1,irep)
               matwk = 0.0d0

               Do im1=istart, iend
                  zsum = 0.0d0
                  Do im2=istart, iend
                     Do iopr=1, nsym
                        z1 = matc(im2,im1,iopr)
                        Do is1=1, ndim_spinor
                           fac = rep_dblegrp(iopr,is1,irep)
                           if ( is1==2 ) z1 = -z1
                           zsum(im2) = zsum(im2) +z1 *fac
                        End Do
                     End Do
                  End Do
                  Do im2=istart, iend
                     norm = zsum(im2) *dble(nd) /dble(nsym*2)
                     if ( abs(norm) > 1.0D-6 ) matwk(im2,im1) = norm
                  ENd Do
               End Do

               call count_basis( immax*ndim_spinor, matwk, num_basis, basis_vec, nadd )

               Do im1=1, nadd
                  orb_irrep(ia,il1,num_basis+im1) = pg_name // ' ' &
                       &                           // symbol_irrep_dblegrp(irep)
               End Do
               num_basis = num_basis +nadd
            End Do
         End Do
         if ( num_basis /= immax*ndim_spinor ) call phase_error_with_msg(nfout,"num_basis /= immax*ndim_spinor"&
                                                                        ,__LINE__,__FILE__)

         allocate( porb_rot_tmp(immax*ndim_spinor,immax*ndim_spinor) )
         porb_rot_tmp(1:immax*ndim_spinor,1:immax*ndim_spinor) &
              &  = porb_rot_matrix_cmplx(ia,il1,1:immax*ndim_spinor,1:immax*ndim_spinor)

         Do im1=1, immax*ndim_spinor
            Do im2=1, immax*ndim_spinor
               ztmp = 0.0d0
               Do im3=1, immax*ndim_spinor
                  ztmp = ztmp +basis_vec(im3,im1) *porb_rot_tmp(im2,im3)
               End Do
               porb_rot_matrix_cmplx(ia,il1,im2,im1) = ztmp
            End Do
         End Do

         deallocate( zsum );      deallocate( matc );  deallocate( work )
         deallocate( basis_vec );  deallocate( porb_rot_tmp )
         deallocate( matwk )
      End Do
      deallocate(crotylm);     deallocate(iylm);     deallocate(nylm)

    end subroutine get_orb_irrep_dblegrp

    subroutine count_basis( ndim, mat, num_basis, basis_vec, nadd )
      integer, intent(in) :: ndim, num_basis
      real(kind=DP), intent(in) :: mat(ndim,ndim)
      integer, intent(out) :: nadd
      real(kind=DP), intent(inout) ::  basis_vec(ndim,ndim)

      integer :: i, lwork, info
      real(kind=DP), allocatable :: rwork(:), eigval(:)

      lwork = 8 *ndim
      allocate( eigval(ndim) );    allocate( rwork(lwork) );
      call dsyev( 'V', 'U', ndim, mat, ndim, eigval, rwork, lwork, info )
      if (info /= 0) then
         write(nfout,*) 'dsyev: info=',info
      end if
      !
      nadd = 0
      Do i=1, ndim
         if ( abs( eigval(i) -1.0 ) < 1.0D-4 ) then
            nadd = nadd +1
            basis_vec(:,num_basis+nadd) = mat(:,i)
         endif
      End do
      deallocate( rwork );  deallocate( eigval )
    end subroutine count_basis

    subroutine calc_local_point_group( max_nn, natom_nn, sym_axis, BondVec, Jtype, &
         &                             nsym0, nsym1, pg_name, rot_space, rot_spinor )
      integer, intent(in) :: max_nn, natom_nn, nsym0
      integer, intent(in) :: Jtype(natom_nn)
      real(kind=DP), intent(in) :: sym_axis(3,3)
      real(kind=DP), intent(in) :: BondVec(max_nn,3)
      integer, intent(out) :: nsym1
      character(len=3), intent(out) :: pg_name
      real(kind=DP), intent(in) :: rot_space(3,3,nsym0)
      complex(kind=CMPLDP), optional, intent(in) :: rot_spinor(2,2,nsym0)

      character(len=5) :: pg_name_i
      integer :: ja, j, k, l, iopr, found, num
      real(kind=DP) :: vec(3), c1, dist
      real(kind=DP), allocatable :: BondVec_rotated(:,:)
      integer, allocatable :: Flag(:), opr_no(:)

      allocate( BondVec_rotated(natom_nn,3) )
      BondVec_rotated= 0.0d0

      Do ja=1, natom_nn
         Do j=1, 3
            c1 = 0.0d0
            Do k=1, 3
               c1 = c1 +sym_axis(k,j) *BondVec(ja,k)
            End Do
            BondVec_rotated(ja,j) = c1
         End Do
      ENd Do
      !
      if ( ipriorb_rot > 2 ) then
         write(nfout,'(A,I5)') "** BondVec of atom :", ia
         Do ja=1, natom_nn
            write(nfout,'(I5,3F20.10)') ja, BondVec(ja,1:3)
         End Do
         write(nfout,'(A)') "BondVec rotated"
         Do ja=1, natom_nn
            write(nfout,'(I5,3F20.10)') ja, BondVec_rotated(ja,1:3)
         End Do
      endif
      !
      allocate( Flag(nsym0) );    Flag = 1
      allocate( opr_no(nsym0) );  opr_no = 0

      Do iopr=1, nsym0
         j_loop: Do j=1, natom_nn
            vec = 0.0d0
            Do k=1, 3
               Do l=1, 3
                  vec(k) = vec(k) +rot_space(k,l,iopr) *BondVec_rotated(j,l)
               End Do
            ENd Do

            found = 0
            k_loop: Do k=1, natom_nn
               dist = sqrt( ( vec(1) -BondVec_rotated(k,1) )**2 &
                    &      +( vec(2) -BondVec_rotated(k,2) )**2 &
                    &      +( vec(3) -BondVec_rotated(k,3) )**2 )
               if ( dist < 1.0E-3 ) then
#if 0
                  if ( Jtype(j) == Jtype(k) ) then
#else
                  if ( Jmass(j) == Jmass(k) ) then
#endif
                     found = 1;   exit k_loop
                  endif
               endif
            End Do k_loop

            if ( found == 0 ) then
#if 0
               if ( natom_nn == 1 .and. Jtype(1) == ityp(ia) ) then
#else
               if ( natom_nn == 1 .and. Jmass(1) == ionic_mass(ia) ) then
#endif
                  c1 = sqrt( (vec(1) +BondVec_rotated(1,1))**2 &
                       &    +(vec(2) +BondVec_rotated(1,2))**2 &
                       &    +(vec(3) +BondVec_rotated(1,3))**2 )
                  if ( c1 < 1.0E-3 ) then
                     found = 1
                  endif
               endif
            endif
100         continue
            if ( found == 0 ) then
               Flag(iopr) = 0
            endif
         ENd Do j_loop
      End Do
       !
      num = 0
      Do iopr=1, nsym0
         if ( Flag(iopr) == 1 ) then
            num = num +1
            opr_no(num) = iopr
            rot_space_wk(:,:,num) = rot_space(:,:,iopr)
            if ( present(rot_spinor) ) rot_spinor_wk(:,:,num) = rot_spinor(:,:,iopr)
         endif
      End Do
      nsym1 = num

      call get_point_group_name( nfout, system, num, opr_no, pg_name, pg_name_i )
      if ( natom_nn == 1 ) then
         if ( pg_name == "D4h" ) pg_name = "Dih"
         if ( pg_name == "C4v" ) pg_name = "Civ"
      endif

      deallocate( BondVec_rotated )
      deallocate( Flag )

    end subroutine calc_local_point_group

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

    subroutine set_symmetry_axis( max_nn, natom_nn, BondVec, axis )
      integer, intent(in) :: max_nn, natom_nn
      real(kind=DP), intent(in) :: BondVec(max_nn,3)

      integer, parameter :: ndim = 3
      real(kind=DP), intent(out) :: axis(ndim,ndim)
      integer :: i, lwork,info, j, num_degene

      real(kind=DP) :: abstol,nfound, ctmp
      real(kind=DP) :: mass, c1, c2, c3, cos_th
      real(kind=DP) :: mat(ndim,ndim)
      real(kind=DP) :: eigval(ndim), eigvec(ndim,ndim)

      integer, allocatable :: iwork(:),ifail(:)
      real(kind=DP), allocatable :: rwork(:)
      real(kind=DP), external :: dlamch
      logical :: flag

      mat = 0.0d0
      Do i=1, natom_nn
         mass = 1.0d0
         mat(1,1) = mat(1,1) +mass *( BondVec(i,2)**2 +BondVec(i,3)**2 )
         mat(2,2) = mat(2,2) +mass *( BondVec(i,3)**2 +BondVec(i,1)**2 )
         mat(3,3) = mat(3,3) +mass *( BondVec(i,1)**2 +BondVec(i,2)**2 )
         mat(1,2) = mat(1,2) -mass *BondVec(i,1) *BondVec(i,2)
         mat(1,3) = mat(1,3) -mass *BondVec(i,1) *BondVec(i,3)
         mat(2,3) = mat(2,3) -mass *BondVec(i,2) *BondVec(i,3)
      End Do
      mat(2,1) = mat(1,2);      mat(3,1) = mat(1,3);      mat(3,2) = mat(2,3)

      lwork = 8 *ndim
      allocate( rwork(lwork) );
      call dsyev( 'V', 'U', ndim, mat, ndim, eigval, rwork, lwork, info )
      !
      if (info /= 0) then
         write(nfout,*) 'dsyev: info=',info
      end if
      deallocate( rwork );
      !
      num_degene = 0
      Do i=1, ndim-1
         Do j=i+1, ndim
            if ( abs(eigval(i) -eigval(j)) < 1.0D-3 ) then
               num_degene = num_degene +1
            endif
         End do
      ENd Do

      if ( ipriorb_rot > 2 ) then
         write(nfout,*) "** number of degeneracy : ", num_degene
         Do i=1, 3
            write(nfout,'(I5,3F20.10)') i, eigval(i)
         End Do
         write(nfout,*) "** eigen_vectors "
         Do i=1, 3
            write(nfout,'(A,I5,A,3F20.10)') "vec ",i, ":", mat(1:3,i)
         End Do
      endif

      select case (num_degene)
      case (0)
         axis(1:3,1:3) = mat(1:3,1:3)
      case (1)
         if ( abs(eigval(1) -eigval(2)) < 1.0D-3 ) then
            axis(1:3,3) = mat(1:3,3)
         else if ( abs(eigval(2) -eigval(3)) < 1.0D-3 ) then
            axis(1:3,3) = mat(1:3,1)
         endif

         flag = .false.
         loop1: Do i=1, natom_nn
            c2 = BondVec(i,1)*axis(1,3) +BondVec(i,2)*axis(2,3) &
                 & +BondVec(i,3)*axis(3,3)
            c3 = sqrt( BondVec(i,1)**2 +BondVec(i,2)**2 +BondVec(i,3)**2 )
            cos_th = c2 /c3
            if ( abs(cos_th) > sqrt(3.)/2. ) cycle
            axis(1:3,1) = BondVec(i,1:3) /c3 -cos_th *axis(1:3,3)

            flag = .true.
            exit loop1
         End Do loop1
         if ( .not. flag ) then
            axis(1,1) = axis(1,3) *2
            axis(2,1) = axis(2,3) *5
            axis(3,1) = -axis(3,3) *7
            c2 = axis(1,1) *axis(1,3) +axis(2,1) *axis(2,3) +axis(3,1) *axis(3,3)
            c3 = sqrt( axis(1,1)**2 +axis(2,1)**2 +axis(3,1)**2 )
            cos_th = c2 /c3
            axis(1:3,1) = axis(1:3,1) /c3 -cos_th *axis(1:3,3)
            flag = .true.
         endif
         axis(1,2) = axis(2,3)*axis(3,1) -axis(3,1)*axis(2,3)
         axis(2,2) = axis(3,3)*axis(1,1) -axis(1,1)*axis(3,3)
         axis(3,2) = axis(1,3)*axis(2,1) -axis(2,1)*axis(1,3)

      case (3)
         c1 = sqrt( BondVec(1,1)**2 +BondVec(1,2)**2 +BondVec(1,3)**2 )
         axis(1:3,1) = BondVec(1,1:3) /c1

         loop2: Do i=2, natom_nn
            c2 = BondVec(i,1)*BondVec(1,1) +BondVec(i,2)*BondVec(1,2) &
                 & +BondVec(i,3)*BondVec(1,3)
            c3 = sqrt( BondVec(i,1)**2 +BondVec(i,2)**2 +BondVec(i,3)**2 )
            cos_th = c2 /c1 /c3
            if ( cos_th < -sqrt(3.)/2. ) cycle
            axis(1:3,2) = BondVec(i,1:3) /c3 -cos_th *axis(1:3,1)
            exit loop2
         End Do loop2
         axis(1,3) = axis(2,1)*axis(3,2) -axis(3,1)*axis(2,2)
         axis(2,3) = axis(3,1)*axis(1,2) -axis(1,1)*axis(3,2)
         axis(3,3) = axis(1,1)*axis(2,2) -axis(2,1)*axis(1,2)
      end select

    end subroutine set_symmetry_axis

    subroutine calc_dist_min( ia, dist_min )
      integer, intent(in) :: ia
      real(kind=DP), intent(out) :: dist_min

      integer :: ja, nx, ny, nz
      real(kind=DP) :: vec(3), dist

      dist_min = 1.0D2
      Loop_ja1: Do ja=1, natm
         Do nx=-1, 1
            Do ny=-1, 1
               Do nz=-1, 1
                  vec(1) = cps(ja,1) -cps(ia,1) &
                       &   +nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                  vec(2) = cps(ja,2) -cps(ia,2) &
                       &   +nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                  vec(3) = cps(ja,3) -cps(ia,3) &
                        &   +nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                  dist = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )
                  if ( dist > 0.01 ) dist_min = min( dist, dist_min )
               ENd Do
            End Do
         End Do
      End Do Loop_ja1
    end subroutine calc_dist_min

    subroutine set_atoms_nn( ia, rcut, max_nn, natom_nn, BondVec, Jtype )
      integer, intent(in) :: max_nn, ia
      real(kind=DP), intent(in) :: rcut
      integer, intent(out) :: natom_nn,  Jtype(max_nn)
      real(kind=DP), intent(out) :: BondVec(max_nn,3)

      integer :: ja, nx, ny, nz
      real(kind=DP) :: vec(3), dist

      natom_nn = 0;   BondVec = 0.0d0;  Jtype = 0

      Loop_ja3: Do ja=1, natm
         Do nx=-1, 1
            Do ny=-1, 1
               Do nz=-1, 1
                  vec(1) = cps(ja,1) -cps(ia,1) &
                       &   +nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                  vec(2) = cps(ja,2) -cps(ia,2) &
                       &   +nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                  vec(3) = cps(ja,3) -cps(ia,3) &
                       &   +nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                  dist = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )

                   if ( dist > 0.01 .and. dist < rcut ) then
                      natom_nn = natom_nn +1
                      BondVec(natom_nn,1:3) = vec(1:3)
                      Jtype(natom_nn) = ityp(ja)
                      Jmass(natom_nn) = ionic_mass(ja)
                   endif
                End Do
             ENd Do
          ENd Do
       End Do Loop_ja3
     end subroutine set_atoms_nn

    subroutine set_porb_rot_matrix_real0( ia )
      integer, intent(in) :: ia

      integer :: il1, val_l, im1, im2, im3, im4
      real(kind=DP) :: c1
      complex(kind=CMPLDP) :: z1, z2, zsum
      complex(kind=CMPLDP), allocatable :: mat(:,:)

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
               porb_rot_matrix_real(ia,il1,im1,im2) = zsum

            End Do
         End Do
         deallocate( mat )
      End Do
    end subroutine set_porb_rot_matrix_real0

    subroutine set_porb_rot_matrix_cmplx0( ia )
      integer, intent(in) :: ia

      integer :: il1, val_l, immax
      integer :: im1, im2, im3, im4, is1, is2, istmp
      real(kind=DP) :: c1
      complex(kind=CMPLDP) :: z1, z2, zsum

      real(kind=DP), allocatable :: rwk_rotmat(:,:)
      complex(kind=CMPLDP), allocatable :: mat(:,:)

      integer :: itmp2, itmp3

      Do il1=1, lmax
         val_l = il1 -1
         immax = 2 *val_l +1
         allocate( mat(-val_l:val_l,-val_l:val_l) );  mat = 0.0d0
         allocate( rwk_rotmat(immax,immax) );  rwk_rotmat = 0.0d0

         Do im1=-val_l ,val_l
            Do im2=-val_l, val_l
               call calc_small_Wigner_function( dble(val_l), dble(im1), dble(im2), &
                    &                           -angle(2), c1 )
               z1 = exp( zi *im1 *angle(1) )
               z2 = exp( zi *im2 *angle(3) )
               mat(im1,im2) = z1 *c1 *z2
            End do
         End do

! ----------------------------------------------------------------
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
               rwk_rotmat(im1,im2) = zsum
            End Do
         End Do
! ----------------------------------------------------------------

         Do im1=1, immax*ndim_spinor
            Do im2=1, immax
               Do is1=1, ndim_spinor
                  zsum = 0.0d0
                  Do im3=1, immax
                     select case (val_l)
                     case (0)
                        zsum = zsum +rwk_rotmat(im3,im2) &
                             &            *JCoeff_L0(immax*(is1-1)+im3,im1)
                     case (1)
                        zsum = zsum +rwk_rotmat(im3,im2) &
                             &            *JCoeff_L1(immax*(is1-1)+im3,im1)
                     case (2)
                        zsum = zsum +rwk_rotmat(im3,im2) &
                             &            *JCoeff_L2(immax*(is1-1)+im3,im1)
                     case (3)
                        zsum = zsum +rwk_rotmat(im3,im2) &
                             &            *JCoeff_L3(immax*(is1-1)+im3,im1)
                     end select
                  End Do
                  porb_rot_matrix_cmplx(ia,il1,immax*(is1-1)+im2,im1) = zsum
               End Do
            End do
         End Do

         deallocate( mat ); deallocate( rwk_rotmat )
      End Do
    end subroutine set_porb_rot_matrix_cmplx0

  end subroutine m_OP_calc_symm

  subroutine set_Jcoeff_L0_L3
    integer :: il1, my_l, immax, num
    integer :: ma, m1, m2, im1
    real(kind=DP) :: val_j, val_mj, c1, c2

    if ( allocated( Jcoeff_L0 ) ) deallocate( Jcoeff_L0 )
    if ( allocated( Jcoeff_L1 ) ) deallocate( Jcoeff_L1 )
    if ( allocated( Jcoeff_L2 ) ) deallocate( Jcoeff_L2 )
    if ( allocated( Jcoeff_L3 ) ) deallocate( Jcoeff_L3 )
    allocate( Jcoeff_L0(1*ndim_spinor,1*ndim_spinor) ); JCoeff_L0 = 0.0d0
    allocate( Jcoeff_L1(3*ndim_spinor,3*ndim_spinor) ); JCoeff_L1 = 0.0d0
    allocate( Jcoeff_L2(5*ndim_spinor,5*ndim_spinor) ); JCoeff_L2 = 0.0d0
    allocate( Jcoeff_L3(7*ndim_spinor,7*ndim_spinor) ); JCoeff_L3 = 0.0d0

    Do il1=1, nloc
       immax = 2*il1 -1;     my_l = il1 -1

       num = 0
! j=l-1/2
       val_j = my_l -0.5d0
       if ( val_j > 0 ) then
          Do ma=-my_l+1, my_l
             val_mj = ma -0.5d0
             m1 = nint(val_mj -0.5d0);   m2 = nint(val_mj +0.5d0)

             c1 = dble( my_l -m1 )/( 2.0d0 *my_l +1.0d0 );  c1 = sqrt(c1)
             c2 = dble( my_l +m2 )/( 2.0d0 *my_l +1.0d0 );  c2 = -sqrt(c2)

             num = num +1
             Do im1=1, immax
                if ( my_l == 1 ) then
                   Jcoeff_L1(      im1,num) = c1 *MatU_ylm_RC_L1(im1,m1)
                   Jcoeff_L1(immax+im1,num) = c2 *MatU_ylm_RC_L1(im1,m2)
                else if ( my_l == 2 ) then
                   JCoeff_L2(      im1,num) = c1 *MatU_ylm_RC_L2(im1,m1)
                   JCoeff_L2(immax+im1,num) = c2 *MatU_ylm_RC_L2(im1,m2)
                else if ( my_l == 3 ) then
                   Jcoeff_L3(      im1,num) = c1 *MatU_ylm_RC_L3(im1,m1)
                   Jcoeff_L3(immax+im1,num) = c2 *MatU_ylm_RC_L3(im1,m2)
                endif
             End Do
          end Do
       endif
! j=l+1/2
       val_j = my_l +0.5d0
       Do ma=-my_l-1, my_l
          val_mj = ma +0.5d0
          m1 = nint(val_mj -0.5d0);  m2 = nint(val_mj +0.5d0)

          c1 = dble( my_l+m2 ) /( 2.0d0 *my_l +1.0d0 );     c1 = sqrt(c1)
          c2 = dble( my_l-m1 ) /( 2.0d0 *my_l +1.0d0 );     c2 = sqrt(c2)

          num = num +1
          if ( m1 >= -my_l .and. m1 <= my_l ) then
             Do im1=1, immax
                if ( my_l == 0 ) then
                   JCoeff_L0(im1,num) = c1 *MatU_ylm_RC_L0(im1,m1)
                else if ( my_l == 1 ) then
                   JCoeff_L1(im1,num) = c1 *MatU_ylm_RC_L1(im1,m1)
                else if ( my_l == 2 ) then
                   JCoeff_L2(im1,num) = c1 *MatU_ylm_RC_L2(im1,m1)
                else if ( my_l == 3 ) then
                   JCoeff_L3(im1,num) = c1 *MatU_ylm_RC_L3(im1,m1)
                endif
             End Do
          end if
          if ( m2 >= -my_l .and. m2 <= my_l ) then
             Do im1=1, immax
                if ( my_l == 0 ) then
                   JCoeff_L0(immax+im1,num) = c2 *MatU_ylm_RC_L0(im1,m2)
                else if ( my_l == 1 ) then
                   JCoeff_L1(immax+im1,num) = c2 *MatU_ylm_RC_L1(im1,m2)
                else if ( my_l == 2 ) then
                   JCoeff_L2(immax+im1,num) = c2 *MatU_ylm_RC_L2(im1,m2)
                else if ( my_l == 3 ) then
                   JCoeff_L3(immax+im1,num) = c2 *MatU_ylm_RC_L3(im1,m2)
                endif
             End Do
          end if
       end Do
    End Do
  end subroutine set_Jcoeff_L0_L3

end module m_OP_rotation
