!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  PROGRAM: TDLRMAIN
!
!  AUTHOR(S): K. Tagami et al   Aug. 1 2011
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================

! ======================================================================
!  This is a collection of subroutines called in the main program.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

!--------------------------------------------------------
!!
!!!           Copy Wfns etc . 
!!
!---------------------------------------------------------
subroutine Copy_Wfn_etc_To_WorkArray( nq )
  use m_LinearResponse_tools,  only : Copy_Zaj_To, Copy_MapZ_To, &
       &                              Copy_NBase_To, Copy_IGf_To, &
       &                              Copy_Iba_To, Copy_Occups_To, &
       &                              Copy_HardPart_To, Copy_EigenVals_To
  use m_LinearResponse_tools,  only : wfn_k, wfn_kmq, map_z_k, map_z_kmq, &
       &                              nbase_k, nbase_kmq, igf_k, igf_kmq, &
       &                              iba_k, iba_kmq, occup_k, occup_kmq, &
       &                              fsval_k, fsval_kmq, &
       &                              eko_k, eko_kmq
  use m_LinearResponse_Control,  only : nrd_efermi, sw_LongWaveLimit
  use m_Const_Parameters,       only : ON, OFF

  use m_PseudoPotential,        only : modnrm
  use m_Const_Parameters,       only : EXECUT  
  use m_Timing,                     only : tstatc0_begin, tstatc0_end

  implicit none
! ----------------------------------------------------
  integer, intent(in)         ::    nq
  integer :: id_sname = -1
! ---------------------------------- start ------------
  call tstatc0_begin('Copy_Wfn_etc_To_WorkArray ', id_sname)
!
  if ( sw_LongWaveLimit == OFF ) then
     if ( nq == 0 ) then
        call Copy_Zaj_To( wfn_k );     call Copy_MapZ_To( map_z_k )
        call Copy_NBase_To( nbase_k ); call Copy_IGf_To( igf_k )
        call Copy_Iba_To( iba_k )
        if ( nrd_efermi==1 )  call Copy_Occups_To( occup_k )
        if ( modnrm==EXECUT ) call Copy_HardPart_To( fsval_k )
     else
        call Copy_Zaj_To( wfn_kmq );     call Copy_MapZ_To( map_z_kmq )
        call Copy_NBase_To( nbase_kmq ); call Copy_IGf_To( igf_kmq )
        call Copy_Iba_To( iba_kmq )
        if ( nrd_efermi==1 ) call Copy_Occups_To( occup_kmq )
        if ( modnrm==EXECUT ) call Copy_HardPart_To( fsval_kmq )
     endif
  else
     if ( nq == 0 ) then
        call Copy_Zaj_To( wfn_k );     call Copy_MapZ_To( map_z_k )
        call Copy_NBase_To( nbase_k ); call Copy_IGf_To( igf_k )
        call Copy_Iba_To( iba_k )
        if ( nrd_efermi==1 ) then
           call Copy_EigenVals_To( eko_k, occup_k )
        else
           call Copy_EigenVals_To( eko_k )
        endif
        if ( modnrm==EXECUT ) call Copy_HardPart_To( fsval_k )
     endif
  endif
  
  call tstatc0_end(id_sname)
  
end subroutine Copy_Wfn_etc_To_WorkArray

!--------------------------------------------------------
!!
!!!   Alloc and Dealloc subroutines
!!
!---------------------------------------------------------
subroutine Alloc_WorkArray_LR1
  use m_LinearResponse_tools,   only : m_LR_alloc_Array_WFns, &
       &                               m_LR_alloc_Array_MapFns,&
       &                               m_LR_alloc_Array_EigenVals, &
       &                               m_LR_alloc_Array_Occups, &
       &                               m_LR_alloc_Array_CorePart
  use m_LinearResponse_Control, only : nrd_efermi
  use m_PseudoPotential,        only : modnrm
  use m_Const_Parameters,       only : EXECUT

  implicit none

  call m_LR_alloc_Array_WFns
  call m_LR_alloc_Array_MapFns
  call m_LR_alloc_Array_EigenVals
  if ( nrd_efermi==1 ) call m_LR_alloc_Array_Occups
  if ( modnrm==EXECUT ) call m_LR_alloc_Array_CorePart
end subroutine Alloc_WorkArray_LR1

subroutine Alloc_WorkArray_LR2
  use m_LinearResponse_Density,   only : m_LR_alloc_Array_Occups_ek, &
       &                               m_LR_alloc_Array_Densities

  implicit none

  call m_LR_alloc_Array_Occups_ek
  call m_LR_alloc_Array_Densities
end subroutine Alloc_WorkArray_LR2

subroutine Free_WorkArray_LR1
  use m_LinearResponse_tools,   only : m_LR_dealloc_Array_WFns, &
       &                               m_LR_dealloc_Array_MapFns,&
       &                               m_LR_dealloc_Array_EigenVals, &
       &                               m_LR_dealloc_Array_Occups, &
       &                               m_LR_dealloc_Array_CorePart
  use m_LinearResponse_Control, only : nrd_efermi
  use m_PseudoPotential,        only : modnrm
  use m_Const_Parameters,       only : EXECUT

  implicit none

  call m_LR_dealloc_Array_WFns
  call m_LR_dealloc_Array_MapFns
  call m_LR_dealloc_Array_EigenVals
  if ( nrd_efermi==1 ) call m_LR_dealloc_Array_Occups
  if ( modnrm==EXECUT ) call m_LR_dealloc_Array_CorePart
  
end subroutine Free_WorkArray_LR1

subroutine Free_WorkArray_LR2
  use m_LinearResponse_Density,   only : m_LR_dealloc_Array_Occups_ek, &
       &                               m_LR_dealloc_Array_Densities

  implicit none
  call m_LR_dealloc_Array_Occups_ek
  call m_LR_dealloc_Array_Densities
end subroutine Free_WorkArray_LR2

!--------------------------------------------------------
!!
!!!   Generate and Delete Q point
!!
!---------------------------------------------------------
subroutine Gen_q_points
  use m_Control_Parameters,       only : Num_q_Points
  use m_LinearResponse_Control,   only : sw_tddft, vqxyz
  use m_LinearResponse_Qpt,       only : q_points_in_BUCS_CARTS

  if ( sw_tddft==1 ) then
     Num_q_Points = 1;
     allocate( vqxyz( Num_q_Points, 3, 2 ) )
     Call q_points_in_BUCS_CARTS
  endif
end subroutine Gen_q_points

subroutine Del_q_points
  use m_LinearResponse_Control,   only : vqxyz
  if ( allocated( vqxyz ) ) Deallocate( vqxyz )
end subroutine Del_q_points

!--------------------------------------------------------
!!
!!!     Initialization and Finalization 
!!
!---------------------------------------------------------
subroutine Initialization_LR
  use m_LinearResponse_Qpt,     only : m_LR_alloc_Vecs_q_plus_G, &
       &                               m_LR_set_Vecs_q_plus_G
  use m_LinearResponse_Qpt,      only : m_LR_alloc_Array_FTQVals, &
       &                                m_LR_calc_FTQ_times_Ylm
  use m_PseudoPotential,         only : modnrm
  use m_Const_Parameters,        only : EXECUT

  use m_Files,                   only : nfout
  use m_LinearResponse_Density,  only : set_ppc_data, &
       &                                Alloc_ptrans, &
       &                                nppcorr, check_PP, &
       &                                Alloc_PP_local_type, &
       &                                Alloc_PP_norm_type, &
       &                                set_ppc_data_it

  use m_Const_Parameters,        only : ON
  use m_LinearResponse_Control,  only : nrd_efermi, sw_LongWaveLimit
  use m_LinearResponse_Qpt,  only : ftqval

! 
  use m_LinearResponse_Control,  only : nmax_G_LR
  use m_PlaneWaveBasisSet,     only : kg0

  implicit none
! ---------------------------- start ----------------
  if ( nrd_efermi==1 ) call Restore_FermiEnergy
  call set_UVSOR_mode_TO( ON )
! 
  if ( nmax_G_LR > kg0 ) nmax_G_LR = kg0
!
  if ( modnrm == EXECUT ) then
     call m_LR_alloc_Array_FTQvals
     call m_LR_calc_FTQ_times_Ylm
  endif
! 
  if ( sw_LongWaveLimit==1 ) then
     call set_nppcorr_to( 2 ) 
     call Alloc_ptrans
!     if ( nppcorr>0 ) call set_ppc_data(nfout)
     if ( nppcorr>0 ) call set_ppc_data_it(nfout)
     call Alloc_PP_norm_type
     call Alloc_PP_local_type
     call check_PP(nfout)
  endif
end subroutine Initialization_LR
   
subroutine Finalization_LR
  use m_LinearResponse_Qpt,  only : m_LR_dealloc_Vecs_q_plus_G, &
       &                            m_LR_dealloc_Array_FTQVals, &
       &                            m_LR_dealloc_Qitg
  use m_PseudoPotential,         only : modnrm
  use m_Const_Parameters,        only : EXECUT
  use m_LinearResponse_Density,  only : Dealloc_Ptrans_Nppc_ilocal_etc, &
       &                                dealloc_PP_local_type, &
       &                                dealloc_PP_norm_type
  use m_LinearResponse_Control,  only : sw_LongWaveLimit
  implicit none
! ------------------------ start ----------------
  call m_LR_dealloc_Vecs_q_plus_G
  if ( modnrm == EXECUT ) then
     call m_LR_dealloc_Qitg
     call m_LR_dealloc_Array_FTQvals
  endif
  if ( sw_LongWaveLimit==1 ) then
     call Dealloc_Ptrans_Nppc_ilocal_etc
     call Dealloc_PP_norm_type
     call Dealloc_PP_local_type
  endif
end subroutine Finalization_LR


subroutine Backup_FermiEnergy
  use m_Electronic_Structure,     only : efermi
  use m_LinearResponse_Control,   only : efermi_backup

  implicit none

  efermi_backup = efermi

end subroutine Backup_FermiEnergy

subroutine Restore_FermiEnergy
  use m_Electronic_Structure,     only : efermi
  use m_LinearResponse_Control,   only : efermi_backup

  implicit none

  efermi = efermi_backup

end subroutine Restore_FermiEnergy

subroutine Set_UVSOR_mode_TO( mode )
  use m_Control_Parameters,       only : uvsormode

  implicit none
  integer, intent(in)  :: mode
  uvsormode = mode

end subroutine Set_UVSOR_mode_TO

subroutine set_nppcorr_to( value )
  use m_LinearResponse_Density, only : nppcorr
  implicit none

  integer, intent(in) :: value
  nppcorr = value

end subroutine set_nppcorr_to

!--------------------------------------------------------
!!
!!!    InputData_Analysis for Linear Response
!!
!--------------------------------------------------------
subroutine InputData_Analysis_LR
  
  use m_Files,                    only : m_Files_reopen_nfinp, nfout
  use m_LinearResponse_Control,   only : m_CtrlP_rd_LinearResponse, &
       &                                 efermi_backup, e_low, e_high, e_step, &
       &                                 nstep, nrd_efermi
  use m_Electronic_Structure,     only : efermi

  implicit none

  call m_Files_reopen_nfinp(1)
  call m_CtrlP_rd_LinearResponse(nfout)
  call m_Files_reopen_nfinp(2)
  call Setup_LinearResponse( nfout, e_low, e_high, e_step, nstep, nrd_efermi, efermi )
  call Backup_FermiEnergy

end subroutine InputData_Analysis_LR

subroutine Setup_LinearResponse( nfout, e_low, e_high, e_step, nstep,  nrd_efermi, efermi )
  use m_Control_Parameters,  only : printable, way_ksample
  use m_Const_Parameters,    only : DP, MESH, Hartree
  use m_LinearResponse_Control, only : e, nistep, scissor
  use m_LinearResponse_Control, only : way_BZintegral, LORENTZIAN_B, GAUSSIAN_B, &
       &                               L_TETRAHEDRON, width, tetra_eps, eta

  implicit none

  integer, intent(in)     :: nfout
  integer, intent(inout)  :: nrd_efermi
  integer, intent(out)    :: nstep
  
  real(DP), intent(in)    :: efermi
  real(DP), intent(inout) :: e_low, e_high, e_step
  !
  integer                 :: i     
  real(DP)                :: e_range
! --------- fermi level ------------
  if(printable) write(nfout,'(1x,"!* nrd_efermi = ",i3)') nrd_efermi
  if(nrd_efermi==0) then
     if(printable) write(nfout,'(1x,"!* efermi is calculated")')
  else
     if(printable) write(nfout,'(1x,"!* efermi = ",f10.5," read from the input file")') efermi
  end if
  ! -------- -- energy range --
  call allocate_energy_range_matrix
  ! ------------ BZ Intelration --
  call BZ_integration_option  
  ! ----------- band gap correction --
  call scissor_operator_option
  
contains 
  
  subroutine allocate_energy_range_matrix
    
    integer i
    !
    if ( e_step <= 0.0d0 ) then
       e_step = 0.002d0
       if(printable) write(nfout,'(1x,"!* default photon energy step of ",f10.5," is used")') e_step
    end if
    ! set default e_low and e_high if necessary
    if ( e_low <= 0.0d0 .and. e_high <= 0.0d0 ) then
       e_low=0.0d0;   e_high=2.0d0
       if ( printable ) then
          write(nfout,'(1x,"!* lowest and highest energy is 0.0d0")')
          write(nfout,'(1x,"!* default lowest energy of ",f10.5," is used ")') e_low
          write(nfout,'(1x,"!* default highest energy of ",f10.5," is used")') e_high
       end if
    end if
    ! set e array
    e_range = ( e_high-e_low ) / e_step
    nstep = int( e_range ) + 1
    
    allocate( e(nstep) ); e = 0.0d0
       
    e(1) = e_low
    if (nstep.gt.1 ) then
       ! 2007.12.10
       do i=2, nstep
          e(i)=e(1)+e_step*(i-1)
       end do
       !         do i=2, nstep
       !            e(i)=e(i-1)+e_step
       !         end do
       ! 2007.10.12
    end if
    if (printable) write(nfout,'(1x,"!* energy range = ",f12.5," -",f12.5,1x,"au",3x," step = ",f12.5,1x,"au")') &
    & e_low,e_high, e_step
  end subroutine allocate_energy_range_matrix
  
  subroutine BZ_integration_option
    if (printable) &
         & write(nfout,'(1x,"!* way_BZintegral = ",i3)') way_BZintegral
    
    if ( printable ) then
       write(nfout,*) '!* eta = ', eta
       write(nfout,*) '!* width = ', width
    endif
    
    if (way_BZintegral==L_TETRAHEDRON) then
       ! ------------------------------------------------------------
!       write(*,*) ' Not supported '
!       stop
       call phase_error_with_msg(nfout, 'Not supported',__LINE__,__FILE__)
       ! ------------------------------------------------------------
       ! check mesh type
       if (way_ksample/=MESH) then
          if(printable) then
             write(nfout,'(1x,"!* linear tetrahedron cannot be used for way_ksample = ",i3)') way_ksample
             write(nfout,'(1x,"!* parabolic broadning is used instead ")')
          end if
          way_BZintegral = LORENTZIAN_B
          call set_default_width
       else
          ! check tetra_eps
          if(printable) &
               &  write(nfout,'(1x,"!* Brillouin zone integration method = linear tetrahedron")')
          if(tetra_eps > 0.0d0) then
             if(printable) write(nfout,'(1x,"!* tetra_eps = ",e12.5," hartree")') tetra_eps
          else
             tetra_eps = 1.0d-4
             if(printable) write(nfout,'(1x,"!* tetra_eps =<0.0d0 is read ; tetra_eps =1.0d-4 hartree is set")')
          end if
          ! check nistep
          if(nistep >= 1) then
             if(nistep > nstep) then
                nistep = nstep
                if(printable) then
                   write(nfout,'(1x,"!* nistep = ",i4," is larger than nstep = ",i4)') nistep, nstep
                   write(nfout,'(1x,"!* nistep = ",i4," is set as default.")') nstep
                end if
             else
                if(printable) write(nfout,'(1x,"!* nistep = ",i4)') nistep
             end if
          else
             if(nistep < 1) then
                nistep = nstep
                if(printable) then
                   write(nfout,'(1x,"!* nistep = ",i4," is set as default.")') nstep
                end if
             end if
          end if
       end if
    end if
    
    if (way_BZintegral/=L_TETRAHEDRON) then
       if (printable) then
          if (way_BZintegral==GAUSSIAN_B) then
             write(nfout,'(1x,"!* Brillouin zone integration method = gaussian broadning")')
             if ( width<=0.0d0 ) then
                call set_default_width
             else
                if(printable) &
                     & write(nfout,'(1x,"!* smearing width = ",f10.5)') width
             end if
          end if
          if (way_BZintegral==LORENTZIAN_B) &
               & write(nfout,'(1x,"!* Brillouin zone integration method = Lorentzian broadning")')
       end if
    end if
    
    !!       if(printable) then
    !!          write(nfout,'(1x,"!* spin = ",i3)') spin
    !!          if(spin/=BOTH) then
    !!             if(spin==MAJOR) write(nfout,'(1x,"!* integration for major spin")')
    !!             if(spin==MINOR) write(nfout,'(1x,"!* integration for minor spin")')
    !!          end if
    !!       end if
  end subroutine BZ_integration_option
  
  subroutine set_default_width
    width=0.5d0 / Hartree  ! 0.01837451d0
    if(printable) &
         & write(nfout,'(1x,"!* default smearing width of 0.01837451 Hartree (=0.5eV) is set")')
  end subroutine set_default_width
  
  subroutine scissor_operator_option
    if(scissor/=0.0d0 .and. printable) then
       write(nfout,'(1x,"!* scissor operator = ", f10.5," Hartree ")') scissor
    end if
  end subroutine scissor_operator_option
  
end subroutine Setup_LinearResponse


