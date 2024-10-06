module m_ES_WF_mixing

  use m_Control_Parameters,  only : nspin, kimg, af, g0_wf_precon, amin_wf_precon, &
       &                            ndim_magmom, noncol
  use m_Const_Parameters,  only : BUCS, DP
  use m_Crystal_Structure,   only : rltv
  use m_KPoints,             only : kv3, vkxyz
  use m_PlaneWaveBasisSet,  only : ngabc, iba, nbase, ngabc, kg1
  use m_Parallelization,     only : MPI_CommGroup, mype, myrank_k, map_k, &
       &                            ista_k, iend_k, np_e

  use m_Electronic_Structure,  only : zaj_l, fsr_l, fsi_l
  use m_ES_nonlocal,       only : m_ES_betar_dot_WFs_4_each_k
  use m_Files,             only : nfout

  use m_Ionic_System,      only : natm
  use m_Charge_Density,   only : hsr, hsro, m_CD_hardpart_hsr, m_CD_hardpart_hsr_noncl
  use m_PseudoPotential,  only : nlmt, nlmta
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP),private,allocatable,dimension(:,:,:,:):: zaj_prev
                                          !d(kg1,np_e,ista_k:iend_k,kimg)

  real(kind=DP), allocatable :: dhsr_precon(:,:,:,:)

contains

  subroutine m_ES_WF_alloc_zaj_prev
    allocate( zaj_prev(kg1,np_e,ista_k:iend_k,kimg) )
  end subroutine m_ES_WF_alloc_zaj_prev

  subroutine m_ES_WF_dealloc_zaj_prev
    deallocate( zaj_prev )
  end subroutine m_ES_WF_dealloc_zaj_prev

  subroutine m_ES_WF_cp_zaj_to_zaj_prev
    integer :: ik, ir

    do ir = 1, kimg
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle                       ! MPI
          zaj_prev(:,:,ik,ir) = zaj_l(:,:,ik,ir)
       end do
    end do
  end subroutine m_ES_WF_cp_zaj_to_zaj_prev

  subroutine m_ES_WF_alloc_dhsr_precon
    allocate( dhsr_precon(natm,nlmt,nlmt,ndim_magmom) )
    dhsr_precon = 0.0d0
  end subroutine m_ES_WF_alloc_dhsr_precon

  subroutine m_ES_WF_dealloc_dhsr_precon
    deallocate( dhsr_precon )
  end subroutine m_ES_WF_dealloc_dhsr_precon

  subroutine m_ES_WF_kerker_mixing( recover_zaj )
    logical, intent(in) :: recover_zaj

    integer :: ik, i, nb
    real(kind=DP) :: ga, gb, gc, gg, q0
    real(kind=DP) :: factor, fac_min
    real(kind=DP) :: ttr(6)

    real(kind=DP), allocatable :: dzaj(:,:,:,:)
    real(kind=DP), allocatable :: zaj_bkup(:,:,:,:)
    real(kind=DP), allocatable :: fsr_bkup(:,:,:), fsi_bkup(:,:,:)

    q0 = g0_wf_precon **2
    fac_min = amin_wf_precon

    allocate( dzaj( kg1, np_e, ista_k:iend_k, kimg ) )
    dzaj = zaj_l -zaj_prev

    if ( recover_zaj ) then
       allocate( zaj_bkup( kg1, np_e, ista_k:iend_k, kimg ) )
       zaj_bkup = zaj_l

       allocate( fsr_bkup(np_e,nlmta,ista_k:iend_k) );
       fsr_bkup = fsr_l

       allocate( fsi_bkup(np_e,nlmta,ista_k:iend_k) );
       fsi_bkup = fsi_l
    endif

    call getttr(rltv,ttr)

    Do ik=1, kv3
       if(map_k(ik) /= myrank_k) cycle

       do i=1, iba(ik)
          nb = nbase(i,ik)
          ga = vkxyz(ik,1,BUCS) + ngabc(nb,1)
          gb = vkxyz(ik,2,BUCS) + ngabc(nb,2)
          gc = vkxyz(ik,3,BUCS) + ngabc(nb,3)
          gg = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               &            + ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
          factor = max( fac_min, gg /(gg+q0) )
          zaj_l(i,:,ik,:) = zaj_prev(i,:,ik,:) + dzaj(i,:,ik,:)*factor
       end do
    End Do

    Do ik=1, kv3
       if(map_k(ik) /= myrank_k) cycle
       call m_ES_betar_dot_WFs_4_each_k(nfout,ik)
    End Do

    if ( noncol ) then
       call m_CD_hardpart_hsr_noncl(nfout,kv3)
    else
       call m_CD_hardpart_hsr(nfout,kv3)
    endif

    deallocate( dzaj )

    if ( recover_zaj ) then
       zaj_l = zaj_bkup
       fsr_l = fsr_bkup;  fsi_l = fsi_bkup
       deallocate( zaj_bkup )
       deallocate( fsr_bkup );  deallocate( fsi_bkup )
    endif

  end subroutine m_ES_WF_kerker_mixing

  subroutine m_ES_WF_precon_hsr
    integer :: ik, i, nb
    real(kind=DP) :: ga, gb, gc, gg, q0
    real(kind=DP) :: factor, fac_min
    real(kind=DP) :: ttr(6)
    real(kind=DP), allocatable :: dzaj(:,:,:,:)
    real(kind=DP), allocatable :: zaj_bkup(:,:,:,:)
    real(kind=DP), allocatable :: hsr_bkup(:,:,:,:)
    real(kind=DP), allocatable :: fsr_bkup(:,:,:), fsi_bkup(:,:,:)

    q0 = g0_wf_precon **2
    fac_min = amin_wf_precon

    allocate( dzaj( kg1, np_e, ista_k:iend_k, kimg ) )
    dzaj = zaj_l -zaj_prev

    allocate( hsr_bkup(natm,nlmt,nlmt,ndim_magmom) )
    hsr_bkup = hsr

    allocate( zaj_bkup( kg1, np_e, ista_k:iend_k, kimg ) )
    zaj_bkup = zaj_l

    allocate( fsr_bkup(np_e,nlmta,ista_k:iend_k) );
    fsr_bkup = fsr_l

    allocate( fsi_bkup(np_e,nlmta,ista_k:iend_k) );
    fsi_bkup = fsi_l

    call getttr(rltv,ttr)

    Do ik=1, kv3
       if(map_k(ik) /= myrank_k) cycle

       do i=1, iba(ik)
          nb = nbase(i,ik)
          ga = vkxyz(ik,1,BUCS) + ngabc(nb,1)
          gb = vkxyz(ik,2,BUCS) + ngabc(nb,2)
          gc = vkxyz(ik,3,BUCS) + ngabc(nb,3)
          gg = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               &            + ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
          factor = max( fac_min, gg /(gg+q0) )
          zaj_l(i,:,ik,:) = zaj_prev(i,:,ik,:) +dzaj(i,:,ik,:)*factor
       end do
    End Do

    Do ik=1, kv3
       if(map_k(ik) /= myrank_k) cycle
       call m_ES_betar_dot_WFs_4_each_k(nfout,ik)
    End Do

    if ( noncol ) then
       call m_CD_hardpart_hsr_noncl(nfout,kv3)
    else
       call m_CD_hardpart_hsr(nfout,kv3)
    endif

    dhsr_precon = hsr -hsro

    zaj_l = zaj_bkup
    hsr = hsr_bkup
    fsr_l = fsr_bkup;  fsi_l = fsi_bkup

    deallocate( hsr_bkup )
    deallocate( fsr_bkup );  deallocate( fsi_bkup )

    deallocate( dzaj );      deallocate( zaj_bkup )

  end subroutine m_ES_WF_precon_hsr

end module m_ES_WF_mixing
