! Copyright (c) 2012, Minoru Otani <minoru.otani@aist.go.jp> 
! 
! Permission is hereby granted, free of charge, to any person 
! obtaining a copy of this software and associated documentation 
! files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, 
! publish, distribute, sublicense, and/or sell copies of the Software, 
! and to permit persons to whom the Software is furnished to do so, 
! subject to the following conditions:
 
! The above copyright notice and this permission notice shall be 
! included in all copies or substantial portions of the Software.
 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
! DEALINGS IN THE SOFTWARE.


Program EsmPack
  Use ESM_VARS
  Implicit None
  Integer :: ig, n1, n2, n3, i
  complex(8),Allocatable :: v_h(:), v_l(:), aux(:)
  Real(8), Allocatable :: force_ew(:,:), force_lc(:,:)

  Character(256) :: inputfile, rhogfile
  Integer :: argc
  Integer :: COMMAND_ARGUMENT_COUNT

  Real(8) :: ewaldr, ewaldg, ewld, ehart

  Real(8), Allocatable :: charge(:)      ! Charge of Atoms
  Real(8) :: alpha
  
  argc = COMMAND_ARGUMENT_COUNT() 
  If( argc < 2 ) Then
     Write(6,*) "Usage : EsmPack.x inputfile rhogfile"
     Stop
  End If

! Setup Variables
  Call GETARG( 1, inputfile )
  Call GETARG( 2, rhogfile  )

  Call SetupEsmVars( inputfile )

  Allocate( charge(nat) )
  charge(:) = upf_zp(:)

  Call ReadRhog( rhogfile )

  Call esm_ggen_2d( )

! Ewald 
  alpha = 2.6d0
  Call Ewald( nat, tau, charge, at, alat, alpha, bg, ewaldr )
  Call esm_ewald( charge, alpha, ewaldg )

  ewld = ewaldg + ewaldr
  
! Hartree
  allocate( v_h(nrxx) )
  Call esm_hartree (rhog_, ehart, v_h )

  write(*,*) "Ewald Energy   : ", ewld
  write(*,*) "Hartree Energy : ", ehart

  allocate( v_l(nrxx) )
  v_l(:) = 0d0
  Call esm_local (v_l)

! Ewald Force
  Allocate( force_ew(3,nat) )
  force_ew(:,:) = 0d0
  alpha = 1d0
  Call esm_force_ew ( alpha, force_ew )
  
! Local Force
  allocate( aux(nrxx) )
  aux(:) = 0d0
  Do ig = 1, ngm
     If( nspin == 1 ) aux(nl(ig)) = rhog_(ig,1)
     If( nspin == 2 ) aux(nl(ig)) = rhog_(ig,1) + rhog_(ig,2)
  End Do
  Allocate( force_lc(3,nat) )
  force_lc(:,:) = 0d0
  Call esm_force_lc ( aux, force_lc )

  write(*,*)
  write(*,*) "Ewald Force : "
  Do i = 1, nat
     write(*,*) force_ew(1:3,i)
  End Do

  write(*,*)
  write(*,*) "Local Force : "
  Do i = 1, nat
     write(*,*) force_lc(1:3,i)
  End Do

  deAllocate( charge, force_ew, force_lc )

  Call esm_printpot( rhog_, v_h, v_l )

  deallocate( v_h, v_l )

  Call DeallocEsmVars

End Program EsmPack

