  subroutine Renew_WF_by_lmSDorlmCG(isolver,precon,dtim)
    integer, intent(in)         :: isolver,precon
    real(kind=DP)               :: dtim
    real(kind=DP), dimension(3) :: etotal
    real(kind=DP), parameter    :: factor = 2
    real(kind=DP)               :: dtim_new = 0.d0, dtim_msdv
    integer                     :: isolver_core, mode

    isolver_core = what_is_the_core_solver(isolver) ! -(in this file)

    mode = ORTHONORMALIZATION

! ---------------- Added by T. Yamasaki, 28 June 2008, 31 Oct. 2008 ---
    if(sw_use_wfred == ON) call m_ESsd_alloc_wfred()
! ---------------------------------<<
    call m_ESsd_copy_zaj_to_zaj_old()

    etotal(1) = m_TE_tell_total_energy()
    dtim_msdv = m_CtrlP_dtim_1Dsearch_now(dtim)
    if(ipri >= 2) write(nfout,'(" !! dtim_msdv = ",f8.4)') dtim_msdv
    if(isolver_core == CG) then
       call m_ESsd_decide_CG_direction(precon)  !  -> betacg
       !    ~~~~~~~~~~~~~~~~~~~~~~~~~~
       call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv)
       !    ~~~~~~~~~~~~~~~~~~~~~~~~~
    else
       call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv)
       !    ~~~~~~~~~~~~~~~~~~~~~~~~~
    end if

    call m_CD_cp_chgq_to_chgqo()
    call ChargeDensity_Construction(1)
    etotal(2) = m_TE_tell_total_energy()

!!$    if(isolver_core == CG .or. isolver_core == eazyCG) mode = NORMALIZATION
    call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,dtim_msdv,factor*dtim_msdv)
	 ! (msdv_grad)
    call ChargeDensity_Construction(1)
    etotal(3) = m_TE_tell_total_energy()
    dtim_new = m_CtrlP_decide_dtim_1Dsearch(nfout,etotal,dtim_msdv,factor)
    mode = ORTHONORMALIZATION
! ------------------- Revised by T. Yamasaki,  03 July 2008, 31 Oct 2008 --->>
    if(sw_use_wfred == ON) then
       call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,dtim_msdv,dtim_new)
    else
       call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,factor*dtim_msdv,dtim_new)
    end if
! ------------------------------------------------------------<<
! ---------------- Added by T. Yamasaki, 28 June 2008, 31 Oct. 2008 ---
    if(sw_use_wfred == ON) call m_ESsd_dealloc_wfred()
! ------------------------------------------------------<<
  end subroutine Renew_WF_by_lmSDorlmCG
