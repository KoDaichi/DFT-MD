!m_container_psicoef::eigenenergyに格納されたエネルギー準位の縮退度を求める
!Copyright (C) 2014-2015 Ryosuke Tomita
! MODIFY SAITO FOR THE CASE ndim_spinor=2 (SO)

module m_degeneracy
#include "defines.F"
  use m_commons, only: nearlyEquals
  use m_container_commons, only: energy_tag, degeneracy, psicoef_t, max_num_deg_band
! SAITO modify **  use m_container_psicoef, only: num_energy_max, UP, DOWN, e_tag
  use m_container_psicoef, only: num_energy_max, UP, DOWN, e_tag, ndim_spinor
  implicit none

  DOUBLE_T, parameter, private:: DEGENERACY_PRECISION = 0.001/27.2116
  integer, private:: num_energy, nspin
!一つのk点における準位の数(縮退したバンドを一つと数える、以下これを縮退点と呼ぶ)
  type(degeneracy), allocatable, target, private::deg(:,:)
!あるk点におけるバンド全体の縮退点の総数(1:nspin)
  integer, private:: num_energy_point(2)
! SAITO 
  integer num_energy_Local,max_num_deg_band_Local
! SAITO
contains
!データのリセット
  subroutine delete_degeneracy
    num_energy_point(1:2) = 0
    deallocate(deg)
  end subroutine delete_degeneracy

!縮退点の数を決定する
  subroutine new_degeneracy(num_energy_t, nspin_t, e_tag_t, is_contains_multirep)
    integer, intent(in):: num_energy_t, nspin_t
!2015/02/16 理由は不明だが、下記のe_tag引数に入っているはずのダウンスピン軌道の固有値がすべて0となってしまう。弥縫策としてcontainer_psicoefのe_tag変数をpublicにして、ここで使用する。


    type(energy_tag), intent(in):: e_tag_t(num_energy_t,nspin_t)
    logical, intent(in):: is_contains_multirep!k点の群が二次元以上の表現をもつかどうか
    integer:: ie, je, is
    
    num_energy = num_energy_t
    nspin = nspin_t
    allocate(deg(num_energy,nspin))
    deg(1:num_energy,1:nspin)%ne = 0
    num_energy_point(1:2) = 0

!一番上の縮退点
#define DEGCUR deg(num_energy_point(is),is)

!k点の群が一次元の表現のみをもつとき、すべての状態が分裂しているものとみなす。
! pairingに対応していないので以下コメントアウト
!    if (ndim_spinor .ne. 2) then  !SAITO modify
!    if(.not.is_contains_multirep) then
!       do is = 1, nspin
!          do ie = 1, num_energy
!縮退点の総数を加算
!             num_energy_point(is) = num_energy_point(is) + 1
!縮退度を1に設定
!             DEGCUR%ne = 1
!縮退点に含まれる準位にieを追加
!             DEGCUR%ie(1) = ie
!          end do
!       end do                    
!       return
!    end if
!    end if                       !SAITO 
! Saito add  ********
    num_energy_Local=num_energy
    max_num_deg_band_Local=max_num_deg_band
    if (ndim_spinor .eq. 2) then
        num_energy_Local=num_energy-1
        max_num_deg_band_Local=4
    endif
!***********************************         
    do is = 1, nspin
       ie = 1
!バンドを1からループ
       do while(ie <= num_energy_Local)
!縮退点の総数を加算
          num_energy_point(is) = num_energy_point(is) + 1
!縮退度を1に設定
          DEGCUR%ne = 1
!縮退点に含まれる準位にieを追加
          DEGCUR%ie(DEGCUR%ne) = ie
          je = ie + 1
!ieより上に存在するバンドのループ
          do while(je <= num_energy)
!ieとjeが縮退しているとき、現在の縮退点の縮退度をインクリメントし、縮退点に含まれる準位にjeを追加する
!ieとjeが縮退していないとき、次のieのループがje+1から始まるようにする
             if((DEGCUR%ne < max_num_deg_band_Local) .and.&
                  isDeg(ie,is,je,is)) then!縮退するとき
                DEGCUR%ne = DEGCUR%ne + 1
                DEGCUR%ie(DEGCUR%ne) = je
                if(je == num_energy) go to 20
                je = je + 1
             else!縮退しないとき
                ie = je
                go to 30
             end if
          end do
          ie = ie + 1
30        continue
       end do
20     continue
    end do
#undef DEGCUR
  contains
!準位(ie,ispin)と(je,jspin)が縮退しているとき.true.、していないとき.false.を返す
    logical function isDeg(ie, ispin, je, jspin)
      integer, intent(in):: ie, ispin, je, jspin
      isDeg = nearlyEquals(e_tag(ie,ispin)%hartree,&
           e_tag(je,jspin)%hartree,&
           DEGENERACY_PRECISION)
    end function isDeg
  end subroutine new_degeneracy

!縮退度の出力
  subroutine print_degeneracy(nfout)
    integer, intent(in):: nfout!出力装置番号
    integer:: ideg, is
    write(nfout, *) '!print_degeneracy'
    do is = 1, nspin
       do ideg = 1, num_energy_point(is)
          write(nfout, fmt = '(2(I0,1X))') ideg, deg(ideg,is)%ne
       end do
    end do
  end subroutine print_degeneracy

!波動関数の構造体に縮退度の情報を追加
  subroutine set_degeneracy(ps_t)
    type(psicoef_t), target, intent(out):: ps_t
    ps_t%num_deg_point = num_energy_point
    ps_t%deg => deg
  end subroutine set_degeneracy
end module m_degeneracy
