!簡約計算をまとめて行うためのモジュール
!Copyright (C) 2014-2015 Ryosuke Tomita

#define READSTOP(func, fileno) \
if(.not.func(fileno)) stop 'error while '//"func"

module m_classification
  use m_commons, only: stdout
  use m_tspace_defines, only: character_table
  use m_tspace_commons, only: print_kpoint_and_kgroup
  use m_container_commons, only: psicoef_t, transd_vector_index, k_state
  use m_container_psicoef, only: read_energy, read_gvector, read_psicoef, get_psicoef, readstate_psicoef, new_psicoef, set_gvector_to_braveis, delete_psicoef, allocate_psicoef, print_psicoef
  use m_transd_vector, only: get_transd_vector_index, new_transd_vector, delete_transd_vector, print_transd_vector
  use m_degeneracy, only: new_degeneracy, set_degeneracy, delete_degeneracy, print_degeneracy
  use m_reduce_band, only: delete_reduce_band, new_reduce_band, print_reduce_band

  type(character_table), target, private:: cht
  type(k_state), target, private:: kst
  type(psicoef_t), target, private:: pst
  type(transd_vector_index), target, private:: tvi
  logical, private:: is_gv_read = .false.
  
contains
!入力されたk点と指標から、電子状態を既約表現に簡約する
  subroutine calcAll(kstt, chtt, fileno)
    integer, intent(in):: fileno!データの出力装置番号
    type(k_state), intent(in):: kstt!k点の情報
    type(character_table), intent(in):: chtt!指標表
    integer:: is, ib
    
    kst = kstt
!指標表の情報を変数に格納
    cht = chtt
!波動関数の準備（空間情報の入力)
    call new_psicoef(kst%kpoint%isdouble)
!エネルギー準位の読み込み
    READSTOP(read_energy, fileno)
!G点の読み込み
    if(.not.is_gv_read) then
       READSTOP(read_gvector, fileno)
       call set_gvector_to_braveis(kst%trp2b)
       is_gv_read = .true.
    end if
!波動関数の準備（変数のアロケートなど）
    call allocate_psicoef
!波動関数の読み込み
    READSTOP(read_psicoef, fileno)
    if(kst%kpoint%isdouble) then
!スピン分極した系のとき、DOWNスピンの状態を読み込む
       READSTOP(readstate_psicoef, fileno)
       READSTOP(read_energy, fileno)
       READSTOP(read_psicoef, fileno)
    end if
    pst = get_psicoef()
! エネルギーの縮退組分け
    call new_degeneracy(pst%num_energy, pst%nspin, pst%e_tag,cht%maxdim > 1)
    call set_degeneracy(pst)
!k点群をもとに、kとGの変換を行う
    call new_transd_vector(kst, pst%gv)
    tvi = get_transd_vector_index()
#ifdef DEBUG
!データのデバッグ
    write(stdout, *) '!Degeneracy'
    call print_degeneracy(stdout)
    write(stdout, *) '!Psicoef'
    call print_psicoef(stdout)
    write(stdout, *) '!transvector'
    call print_transd_vector(stdout)
#endif
!可約表現の計算
    call new_reduce_band(cht, kst, pst, tvi)
  end subroutine calcAll

!データの出力
  subroutine printAll(kstt, is, fno)
    type(k_state), intent(in):: kstt!k点の情報
    integer, intent(in):: fno, is!出力装置番号とスピンインデックス
    write(fno, fmt = '(A)', advance = 'no') '!'
    call print_kpoint_and_kgroup(fno, kstt%kpoint, kstt%kgroup)
!簡約係数の出力
    call print_reduce_band(fno, is)
  end subroutine printAll

!アロケートしたデータのリセット
  subroutine deleteAll
    call delete_reduce_band
    call delete_degeneracy
    call delete_psicoef
    call delete_transd_vector
  end subroutine deleteAll
end module m_classification
