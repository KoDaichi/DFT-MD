!PHASEから入力されたG点の座標、エネルギー固有値、波動関数を読み込み、保存するモジュール
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_container_psicoef
#include "defines.F"
  use m_container_commons, only: energy_tag, psicoef_t, linebuf, parser_tag, readline, readline_until_tag
  implicit none

!m_container_common::readline_until_tag用のタグ
  character(*), parameter, private::&
       gvector_tag = 'GVECTOR',&
       eigen_energy_tag = 'EIGEN_ENERGY',&
       psicoef_tag = 'PSICOEF',&
       magnetic_state_tag = 'MAGNETIC_STATE',&
       spinstate_tag = 'SPINSTATE',&
       numgv_tag = 'Number_of_Gvectors',&
       eigen_energy_num_tag = 'EIGEN_ENERGY NUMBER',&
       kpoint_tag = 'KPOINT'

!エネルギー準位およびGベクトルの数、後の関数で指定される
  integer, private:: num_energy = 0, num_gvector = 0
!エネルギー準位及びGベクトルの最大数
  integer, parameter:: num_energy_max = 4095, num_gvector_max = 65535, UP = 1, DOWN = 2
  integer, parameter, private:: NUM_COEF_COLUMN = 4
!Gベクトルの成分
  integer, private, target:: gvector(3, num_gvector_max)
  type(energy_tag), public, target:: e_tag(num_energy_max,2)
!波動関数の成分psicoef(1:num_gvector, 1:num_energy, 2)、領域は後の関数で確保される、psicoef(:,:,1)はアップ軌道、2はダウン軌道
  COMPLEX_T, allocatable, private, target:: psicoef(:,:,:)
  integer, public:: nspin = 1, af = 0, ikpoint = 1, ispin = 1, ndim_spinor = 1
contains
!波動関数のデータ読み込み
  logical function read_psicoef(nfin)&
       result(isread)
    integer, intent(in):: nfin
    integer:: ie, ig, irow, icol
    DOUBLE_T:: c(2,NUM_COEF_COLUMN)
    DOUBLE_T:: nrm
    isread = .false.
!new_psicoefが実行されていないのでFALSEを返す
    if(.not.allocated(psicoef)) return
!読み込み開始処理、波動関数タグまで読み飛ばし、タグがなければ異常終了
    if(readline_until_tag(nfin, psicoef_tag) < 1) return
!読み込み終了処理、読み込み終了タグもしくはEOFまで読み込む、読み込みエラーの時は異常終了

    ie = 1
    ig = 0
    do while(.true.)
       READSELECT(nfin)
       if(index(linebuf, eigen_energy_num_tag) > 0) then
          read(linebuf, *) ie
          ig = 0
          cycle
       end if
       if(index(linebuf, parser_tag) > 0) exit
       if(ie > num_energy .or. ig >= num_gvector) cycle
       read(linebuf, *) (c(1:2, icol), icol = 1,NUM_COEF_COLUMN)
       do icol = 1, NUM_COEF_COLUMN
          psicoef(ig + icol, ie, ispin) = CMPLX_(c(1, icol), c(2, icol))
       enddo
       ig = ig + NUM_COEF_COLUMN
    end do
    isread = .true.
!規格化
    do ie = 1, num_energy
       nrm = 0d0
       do ig = 1, num_gvector
          nrm = nrm + real(conjg(psicoef(ig, ie, ispin))* psicoef(ig, ie, ispin))
       end do
       nrm = 1d0/sqrt(nrm)
       do ig = 1, num_gvector
          psicoef(ig, ie, ispin) = psicoef(ig, ie, ispin) * nrm
       end do
    end do
  end function read_psicoef

!エネルギー準位および数値を読み込み、準位数num_energyを決定する
  logical function read_energy(nfin) result(isread)
    integer, intent(in):: nfin
    integer:: ie
    DOUBLE_T:: ene

    isread = .false.
!読み込み開始処理、psicoefと同様
    if(readline_until_tag(nfin, eigen_energy_tag) < 1) return
!読み込み終了処理、psicoefと同様
    do while(.true.)
       READSELECT(nfin)
       if(index(linebuf, parser_tag) > 0) exit
       read(linebuf, *) ie, ene
       if(ie > num_energy_max) cycle
       e_tag(ie,ispin)%ie = ie
       e_tag(ie,ispin)%hartree = ene
!最大数更新、num_energyの値を更新
       num_energy = max(ie, num_energy)
    end do

    isread = .true.
  end function read_energy

!磁性状態(スピンの数とantiferroかどうか)を読み込む
  logical function readmagstate_psicoef(nfin) result(is)
    integer, intent(in):: nfin
    is = .false.
!magnetic_state_tagまで読み飛ばし
    if(readline_until_tag(nfin, magnetic_state_tag) < 1 .or. &
         readline(nfin) < 1) return
!-----------------------------Modified by T.A.Ariasoca-----------------------------
    read(linebuf, *) nspin, af, ndim_spinor
!---------------------------------------------------------------------------
    is = .true.
  end function readmagstate_psicoef

!k点のインデックスとスピン状態を読み込む
  logical function readstate_psicoef(nfin) result(is)
    integer :: i
    integer, intent(in):: nfin
    DOUBLE_T:: dummy(3)
    character(len=4):: spinstate
    is = .false.
!kpoint_tagまで読み飛ばし
    if(readline_until_tag(nfin, kpoint_tag) < 1 .or. &
         readline(nfin) < 1) return
!-----------------------------Modified by T.A.Ariasoca-----------------------------
    if(nspin == 2 .and. ndim_spinor == 2) then
       read(linebuf, *) ikpoint
       if(mod(ikpoint,2)==1) then
          ispin = 1
          ikpoint = ikpoint/2 + 1
       else
          ispin = 2
          ikpoint = ikpoint/2
       end if
!----------------------------------------------------------------------------------
    else if(nspin == 2) then
       read(linebuf, *) ikpoint, dummy(1:3), spinstate
       if(spinstate(1:2) == 'UP') then !アップスピンのとき
          ispin = 1
          ikpoint = ikpoint/2 + 1
       else!ダウンスピンのとき
          ispin = 2
          ikpoint = ikpoint/2
       end if
    else!非磁性のとき
       read(linebuf, *) ikpoint
    end if
    is = .true.
  end function readstate_psicoef

!Gベクトルの成分の読み込み
  logical function read_gvector(nfin)
    integer, intent(in):: nfin
    integer:: a, gbuf(3)

    read_gvector = .false.
    if(readline_until_tag(nfin, numgv_tag) < 1 .or. &
         readline(nfin) < 1) return
    read(linebuf, *) num_gvector

!読み込み開始処理、psicoefと同様
    if(readline_until_tag(nfin, gvector_tag) < 1) return
!読み込み終了処理、psicoefと同様
    do while(.true.)
       READSELECT(nfin)
       if(index(linebuf, parser_tag) > 0) exit
       read(linebuf, *) a, gbuf
!領域外処理, psicoefと同様
       if(a > num_gvector) cycle
       gvector(:, a) = gbuf
    end do
    read_gvector = .true.
  end function read_gvector

!Gベクトル(基本逆格子の成分)をブラベ逆格子の成分に変換
  subroutine set_gvector_to_braveis(trp2b)
    DOUBLE_T, intent(in):: trp2b(3,3)
    integer:: ig
    do ig = 1, num_gvector
       gvector(1:3,ig) = nint(matmul(dble(gvector(1:3,ig)), trp2b))
    enddo
  end subroutine set_gvector_to_braveis

!得られたデータを構造体にコピーし、それを返す
  type(psicoef_t) function get_psicoef() result(p)
    p%num_energy = num_energy
    p%gv%num = num_gvector
    p%gv%comp => gvector
    p%nspin = nspin
    p%e_tag => e_tag
    p%comp => psicoef
!--------------------Modified by T.A.Ariasoca---------------------------------    
  ! p%nspinor = ndim_spinor
!-----------------------------------------------------------------------------
  end function get_psicoef

!k点の情報をコピーする(スピンの数およびk点の格子変換行列)
  subroutine new_psicoef(isdouble)
    logical, intent(in):: isdouble
    nspin = 1
!--------------------Modified by T.A.Ariasoca-------------------------------------
    !if(isnoncol) nspin = 2
!--------------------------------------------------------------------------------
    if(isdouble) nspin = 2
  end subroutine new_psicoef

!read_energy, read_gvectorで設定した変数num_energy, num_gvectorをもってpsicoefの領域を確保
  subroutine allocate_psicoef
    if(.not.allocated(psicoef))&
         allocate(psicoef(num_gvector, num_energy, nspin))
  end subroutine allocate_psicoef

!psicoefの領域を開放し、データをリセットする
  subroutine delete_psicoef
    num_energy = 0
!    num_gvector = 0
    if(allocated(psicoef)) deallocate(psicoef)
  end subroutine delete_psicoef

!設定したデータの出力
  subroutine print_psicoef(nfout)
    integer, intent(in):: nfout
    integer:: ie, ig, is

    write(nfout, *) num_energy, num_gvector, nspin

    do ig = 1, num_gvector
       write(nfout, *) ig, gvector(:, ig)
    enddo
    write(nfout, *)
    do is = 1, nspin
       write(nfout, *) 'ISPIN:', is
       do ie = 1, num_energy
          write(nfout, *) 'ENERGY:', e_tag(ie,is)%ie, e_tag(ie,is)%hartree
          do ig = 1, num_gvector
             write(nfout, *) ig, psicoef(ig, e_tag(ie,is)%ie, 1:nspin)
          enddo
          write(nfout, *)
       enddo
    end do
  end subroutine print_psicoef

end module m_container_psicoef
