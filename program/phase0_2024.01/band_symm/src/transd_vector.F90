!kベクトルおよびGベクトルに対称操作を施す処理
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_transd_vector
#include "defines.F"
#define ROTATION(rot, vec) matmul(vec, rot)
!$ use omp_lib
  use m_commons, only: equals_vector
  use m_container_commons, only: k_state, gvector_t, transd_vector_index
  use m_tspace_defines, only: group, num_operator_max

  implicit none
!transd_gv_index(1: kgroup%order): G点を回転操作Rにより変換させたときの新たなG'のインデックス。(*)そのG点がm_container_psicoef::gvectorにあるとき、そのインデックスを、ないとき-1を格納する
!transd_kv_index(1: kgroup%order): k点がBZ境界にあるとき、等価な点G = k - Rkが存在する。そのG点が~(上記*と同じ)~格納する
  type(k_state), pointer, private:: ks
  type(gvector_t), pointer, private:: gv
  integer, private, allocatable, target:: &
       transd_gv_index(:,:), transd_kv_index(:)
contains
!Gベクトルをgv_t%compから検索し、そのインデックスを返す、みつからなかったら-1が返る
  integer function find_gvector(kv, gv_t)
    integer, dimension(3), intent(in):: kv
    type(gvector_t), intent(in):: gv_t
    integer:: a
    find_gvector = -1
    do a = 1, gv_t%num
       if(equals_vector(kv, gv_t%comp(1:3, a))) then
          find_gvector = a
          exit
       endif
    enddo
  end function find_gvector

!指定したk点のtransd_kv_indexを返す
  integer function find_rotk(iopr)
    integer, intent(in):: iopr
    integer:: a, n, itmp(3), jtmp(3)

!任意のk点のインデックスは有理数、つまり共通分母整数Nと分子整数ベクトルG'で表される, k = G'/N
!G' - RG'  = Ng (defined as g, k - Rk = g)
    itmp = ks%kpoint%comp - nint(ROTATION(ks%kgroup%op(iopr)%rot, ks%kpoint%comp))
    find_rotk = -1
    n = ks%kpoint%common_denom
!g == G <=> G' == NG
    do a = 1, gv%num
       jtmp = n * gv%comp(1:3, a)
       if(equals_vector(itmp, jtmp)) then
          find_rotk = a
          exit
       endif
    enddo
  end function find_rotk

  subroutine delete_transd_vector
    if(allocated(transd_kv_index))&
         deallocate(transd_kv_index)
    if(allocated(transd_gv_index))&
         deallocate(transd_gv_index)
  end subroutine delete_transd_vector

!k点群の対称操作を指定する（任意）
  subroutine new_transd_vector(ks_t, gv_t)
    type(k_state), target, intent(in):: ks_t
    type(gvector_t), target, intent(in):: gv_t
    integer:: iopr, ig, itmp(3)
    DOUBLE_T:: rot(3,3)

    ks => ks_t
    gv => gv_t
!k点がBZ境界にあるとき、等価な点の処理が必要
    if(.not.allocated(transd_kv_index) .and. ks%kpoint%comm2bz > 1) then
       allocate(transd_kv_index(ks%kgroup%order))
       do iopr = 1, ks%kgroup%order
          transd_kv_index(iopr) = find_rotk(iopr)
       end do
    endif
    if(.not.allocated(transd_gv_index)) then
       allocate(transd_gv_index(gv%num, ks%kgroup%order))
    else
       return
    end if

!$omp parallel
    do iopr = 1, ks%kgroup%order
       rot = ks%kgroup%op(iopr)%rot
!$omp do private(itmp)
       do ig = 1, gv%num
          itmp = nint(ROTATION(rot,gv%comp(1:3,ig)))! G' = RG
          transd_gv_index(ig, iopr) = find_gvector(itmp, gv)
       end do
!$omp end do
    end do
!$omp end parallel
  end subroutine new_transd_vector

!変換されたベクトルの出力
  subroutine print_transd_vector(nfout)
    integer, intent(in):: nfout
    integer:: iopr, ig
    do iopr = 1, ks%kgroup%order
       write(nfout, *) 'K, ', ks%kgroup%op(iopr)%name
       if(ks%kpoint%comm2bz > 1)&
            write(nfout, *) transd_kv_index(iopr)
       write(nfout, *) 'G, ', ks%kgroup%op(iopr)%name
#ifdef NUM_G_DEBUG
       do ig = 1, NUM_G_DEBUG
#else
       do ig = 1, gv%num
#endif
         write(nfout, *) ig, transd_gv_index(ig, iopr)
       end do
    end do
  end subroutine print_transd_vector

  type(transd_vector_index) function &
       get_transd_vector_index() result(tvi)
  tvi%numop = ks%kgroup%order
  tvi%numg = gv%num
  tvi%gv => transd_gv_index
  tvi%kv => transd_kv_index
  end function get_transd_vector_index

end module m_transd_vector
