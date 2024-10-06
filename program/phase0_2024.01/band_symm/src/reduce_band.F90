!PHASEから読み込んだ平面波波動関数とk群の射影演算子との計算を行うモジュール
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_reduce_band
!$ use omp_lib
!読み込んだデータ及びTSPACEから得られる情報から、波動関数を簡約する
#include "defines.F"
#ifdef NUM_G_DEBUG
#define G_DEBUG ig<NUM_G_DEBUG
#else
#define G_DEBUG .true.
#endif

  use m_commons, only: pi2, stdout, zi
  use m_tspace_defines, only: irrep, character_table
  use m_container_commons, only: k_state, psicoef_t, transd_vector_index
!---------------------Modified by T.A.Ariasoca-----------------------------
!import ndim_spinor from container_psicoef module
  use m_container_psicoef, only: ndim_spinor
!--------------------------------------------------------------------------
  use m_character_table, only: get_trace
  use m_transd_vector, only: find_gvector
  implicit none

!kベクトルの実数表現
  DOUBLE_T, dimension(nxyz_), private:: kv
!現在のk点バンドの可約表現(redrep(縮退エネルギー準位,スピンインデックス))
  type(irrep), allocatable, target, private:: redrep(:,:)
  logical, private:: is_special_k

  type(k_state), pointer, private:: ks
  type(psicoef_t), pointer, private:: ps
  type(transd_vector_index), pointer, private:: tvi
  type(character_table), pointer, private:: chtable
contains
!特定のエネルギー準位、対称操作における可約表現の行列成分をもとめる(<psi^k_ie1,is1|R|psi^k_ie2,is2>)
  function get_reducible_rep(ie1, ie2, is1, is2, iopr)&
       result(sbuf)
    integer, intent(in):: ie1, ie2, is1, is2, iopr

    integer:: ig, new_ig
    integer, dimension(nxyz_):: new_ig_vec, Rk !RG - g
        COMPLEX_T:: cc, sbuf
    DOUBLE_T:: ktr
#define IS_OUTER_BOUNDS(ig) (ig < 1 .or. ig > ps%gv%num)

    sbuf = 0d0
    Rk = 0
    if(is_special_k) then
       new_ig = tvi%kv(iopr)
!領域外処理
       if(.not.IS_OUTER_BOUNDS(new_ig)) then
          Rk = ps%gv%comp(:, new_ig)
       else
          return
       end if
    end if
    ktr = dot_product(- Rk, ks%kgroup%op(iopr)%tr)
#ifdef DEBUG
    write(stdout, *) '#make_psi_op_psi#kp',ks%kpoint%name,'ie1 ', ie1, 'ie2 ', ie2, 'opr ', ks%kgroup%op(iopr)%name, 'iopr', iopr
    write(stdout, *) 'kv ', kv, 'RK ', Rk
    write(stdout, *) 'tr ', ks%kgroup%op(iopr)%tr
#endif

!$omp parallel
!$omp do private(new_ig, new_ig_vec, cc), reduction(+:sbuf)
!全g点に関する処理
    do ig = 1, ps%gv%num
!対称操作により変化したg点のインデックスをコピー
       new_ig = tvi%gv(ig, iopr)
!領域外処理
       if(IS_OUTER_BOUNDS(new_ig)) cycle

#ifdef DEBUG
       if(G_DEBUG) then
          write(stdout, *) 'Transd_gv_index', ig, ps%gv%comp(:, ig)
          write(stdout, *) 'Transd_gv_index', new_ig, ps%gv%comp(:, new_ig)
          write(stdout, *)
       end if
#endif
       new_ig_vec = ps%gv%comp(:, new_ig)! RG

!k点が非共型であり、BZ境界のk点のとき
       if(is_special_k) then
          new_ig = find_gvector(new_ig_vec - Rk, ps%gv)
!領域外処理
          if(IS_OUTER_BOUNDS(new_ig)) cycle
       end if

!並進操作からなる位相cc(k点がBZ境界にあるとき、1でない値を持つ、k点がBZ内部にあるときは常に１)
       cc = CEXP_(-pi2 * (ktr + dot_product(new_ig_vec, ks%kgroup%op(iopr)%tr))*zi)

#define RIGHT ps%comp(ig, ie2, is2)
#define LEFT ps%comp(new_ig, ie1, is1)

#ifdef DEBUG
       if(G_DEBUG) then
          write(stdout, *) 'psi(', ig, ')= ', RIGHT
          write(stdout, *) 'psi(', new_ig, ')= ', LEFT
          write(stdout, *) 'k-g+RG ', kv + new_ig_vec, 'tr ', ks%kgroup%op(iopr)%tr
          write(stdout, *) 'cc ', cc, sbuf
          write(stdout, *)
       endif
#endif
       sbuf = sbuf + cc*conjg(LEFT)*RIGHT
    enddo
#undef RIGHT
#undef LEFT
!$omp end do
!$omp end parallel
#ifdef DEBUG
    write(stdout, *) '<psi|op|psi>= ', sbuf
    write(stdout, *)
#endif
#undef IS_OUTER_BOUNDS
  end function get_reducible_rep

!---------------------Modified by T.A.Ariasoca----------------------------------------
!Subroutine to calculate double group irreducible representation, modified version of get_reducible_rep subroutine
  function get_reducible_rep_noncol(ie1, ie2, is1, iopr)&
       result(sbuf)
    integer, intent(in):: ie1, ie2, is1, iopr

    integer:: ig, new_ig
    integer, dimension(nxyz_):: new_ig_vec, Rk !RG - g
        COMPLEX_T:: cc, sbuf, dg1, dg2
    DOUBLE_T:: ktr
#define IS_OUTER_BOUNDS(ig) (ig < 1 .or. ig > ps%gv%num)

    sbuf = 0d0
    Rk = 0
    if(is_special_k) then
       new_ig = tvi%kv(iopr)
!if special kpoint
       if(.not.IS_OUTER_BOUNDS(new_ig)) then
          Rk = ps%gv%comp(:, new_ig)
       else
          return
       end if
    end if
    ktr = dot_product(- Rk, ks%kgroup%op(iopr)%tr)

!$omp parallel
!$omp do private(new_ig, new_ig_vec, cc), reduction(+:sbuf)
    do ig = 1, ps%gv%num
       new_ig = tvi%gv(ig, iopr)

       if(IS_OUTER_BOUNDS(new_ig)) cycle

       new_ig_vec = ps%gv%comp(:, new_ig)! RG

       if(is_special_k) then
          new_ig = find_gvector(new_ig_vec - Rk, ps%gv)

          if(IS_OUTER_BOUNDS(new_ig)) cycle
       end if

!cc is the translation vector part exp(-i(c_i*G_n)\dot\tau_i)
       cc = CEXP_(-pi2 * (ktr + dot_product(new_ig_vec, ks%kgroup%op(iopr)%tr))*zi)

!dg1 and dg2 correspond to spin rotation
       dg1 = (ks%kgroup%op(iopr)%spinrot(1,1)*ps%comp(ig, ie2, is1))+(ks%kgroup%op(iopr)%spinrot(1,2)*ps%comp(ig, ie2, is1+1))
       dg2 = (ks%kgroup%op(iopr)%spinrot(2,1)*ps%comp(ig, ie2, is1))+(ks%kgroup%op(iopr)%spinrot(2,2)*ps%comp(ig, ie2, is1+1))

!sbuf is the modified version of <psi|op|psi> with double group operation
!bag fixed M.saito
!       sbuf = sbuf + cc*((conjg(ps%comp(new_ig, ie1, is1))*dg1)+(conjg(ps%comp(new_ig, ie1, is1+1)*dg2)))
        sbuf = sbuf + cc*( conjg(ps%comp(new_ig, ie1, is1))*dg1 +conjg(ps%comp(new_ig, ie1, is1+1))*dg2)
    enddo
!$omp end do
!$omp end parallel
#ifdef DEBUG
    write(stdout, *) '<psi|op|psi>= ', sbuf
    write(stdout, *)
#endif
#undef IS_OUTER_BOUNDS
  end function get_reducible_rep_noncol
!---------------------------------------------------------------------------------------------


!get_reducible_repで計算した可約表現の指標から、指定した既約表現imrの係数を求める1/h sum_R x_(R)<psi|R|psi>
  COMPLEX_T function irreducible_coefficient(ideg, is, imr)&
       result(coef)
!エネルギー準位、既約表現インデックス、スピンインデックス
    integer, intent(in):: ideg, imr, is
    integer:: ic
    coef = 0d0

#define CLASS ks%kgroup%cls(ic)
#define COEF(imr) imr%comp(ic)
    do ic = 1, chtable%nclass!ks%kgroup%nclass
       coef = coef + CLASS%order *&
            conjg(COEF(chtable%imrs(imr))) *&
            COEF(redrep(ideg,is))
    end do
    coef = coef / CMPLX_(ks%kgroup%order, 0d0)
#undef CLASS
#undef IRREP
  end function irreducible_coefficient

#define I_ENERGY(id) ps%e_tag(ps%deg(ideg,is)%ie(id),is)
!バンドの可約表現の出力
  subroutine print_reduce_band(nfout,is)
    integer, intent(in):: nfout, is
    integer:: ideg, imr
    do ideg = 1, ps%num_deg_point(is)
       write(nfout, fmt = '(A)', advance = 'no') '!'
!エネルギーのインデックス、縮退点のエネルギー、縮退度
       write(nfout, '(I0,1X,f16.8,1X,I0)') I_ENERGY(1)%ie, I_ENERGY(1)%hartree, ps%deg(ideg,is)%ne
!可約表現の指標及び行列表現の出力
       !call print_cr_wd(ideg, is, nfout)
!可約表現に含まれる、既約表現の係数の出力
       do imr = 1, chtable%numimr
          write(nfout, fmt = '(A,1X)', advance = 'no') trim(chtable%imrs(imr)%name)//'('//trim(chtable%imrs(imr)%name_mulliken)//')'
          write(nfout, 400) irreducible_coefficient(ideg, is, imr), chtable%imrs(imr)%dim
       end do
       write(nfout, *)
    end do
400 format('(',f0.8,',',f0.8,')',1X,I0)
  end subroutine print_reduce_band

  subroutine new_reduce_band(cht, kst, pst, tvit)
    type(k_state), target, intent(in):: kst
    type(psicoef_t), target, intent(in):: pst
    type(character_table), target, intent(in):: cht
    type(transd_vector_index), target, intent(in):: tvit
    integer:: ideg, ic, iopr_ind, is
!k点、波動関数、指標表、kベクトル変化の情報を受け取る(ポインタによるシャローコピー)
    ks => kst
    ps => pst
    chtable => cht
    tvi => tvit
    is_special_k = (ks%kpoint%comm2bz > 1)
    kv(1:3) = dble(ks%kpoint%comp(1:3))/ks%kpoint%common_denom
    if(.not.AllocateAll()) return

!--------------------Modified by T.A.Ariasoca---------------------------------
!Used to get double group irreducible representation, by not looping through nspin
    if(ndim_spinor == 2) then
        is = 1
       do ideg = 1, ps%num_deg_point(is)
          do ic = 1, ks%kgroup%nclass
             call set_cr(ideg, is, ic)
          end do
       end do
    else
!for collinear calculation, loop through all spin
!----------------------------------------------------------------------------
       do is = 1, ps%nspin
!全縮退エネルギー準位におけるループ
          do ideg = 1, ps%num_deg_point(is)
!k点の属する群の全類に関するループ
             do ic = 1, ks%kgroup%nclass
                call set_cr(ideg, is, ic)
             end do
          end do
       end do
    end if

  contains
!指定されたエネルギー準位と対称操作における、可約表現の指標を返し、行列表現を変数に保存する
    subroutine set_cr(ideg, is, ic)
      integer, intent(in):: ideg, ic, is
      integer:: iopr, id, jd!, is, js
      COMPLEX_T:: buf!(2,2)

#define CLS ks%kgroup%cls(ic)
!一つの縮退エネルギー準位に含まれるバンドに対するループ
      redrep(ideg,is)%dim = ps%deg(ideg,is)%ne
      do iopr = 1, CLS%order
         do id = 1, ps%deg(ideg,is)%ne
!同上
            do jd = 1, ps%deg(ideg,is)%ne
!指定された2つのバンドと、対称操作における可約表現の行列表現を返す
!--------------------Modified by T.A.Ariasoca---------------------------------
!Calculate buf for double group representation
               if(ndim_spinor == 2) then
                  buf = get_reducible_rep_noncol(I_ENERGY(id)%ie, I_ENERGY(jd)%ie, is, CLS%ind(iopr))
               else
!-----------------------------------------------------------------------------
                  buf = get_reducible_rep(I_ENERGY(id)%ie, I_ENERGY(jd)%ie, is, is, CLS%ind(iopr))
               end if
               redrep(ideg,is)%wd(id,jd,CLS%ind(iopr)) = buf
            end do
         end do
      end do

      redrep(ideg,is)%comp(ic) = get_trace(&
           redrep(ideg,is)%wd(:,:, CLS%ind(1)), ps%deg(ideg,is)%ne)
#undef CLS
    end subroutine set_cr

    logical function AllocateAll()
      integer:: ndeg
      AllocateAll = .false.
      if(.not.allocated(redrep)) then
         ndeg = 0
         if(ps%num_deg_point(1) > ps%num_deg_point(2)) then
            ndeg = ps%num_deg_point(1)
         else
            ndeg = ps%num_deg_point(2)
         end if
         allocate(redrep(ndeg,ps%nspin))
      end if
      AllocateAll = .true.
    end function AllocateAll
  end subroutine new_reduce_band

  subroutine delete_reduce_band
    if(.not.DeallocateAll()) return
  contains
    logical function DeallocateAll()
      DeallocateAll = .false.
      if(allocated(redrep)) deallocate(redrep)
      DeallocateAll = .true.
    end function DeallocateAll
  end subroutine delete_reduce_band

!指定されたエネルギー準位の可約表現を出力
  subroutine print_cr_wd(ideg, is, nfout)
    integer, intent(in):: ideg, is, nfout
    integer:: ic, iopr_ind, id, jd, dim
!可約表現の次元。バンドの縮退数
    dim = ps%deg(ideg,is)%ne
    write(nfout, *) '##reducible representation ', ks%kpoint%name
!類に関するループ
    do ic = 1, ks%kgroup%nclass
!類の名前及び可約表現の指標の出力
       write(nfout, 401) '##', ks%kgroup%cls(ic)%name, redrep(ideg,is)%comp(ic)
!類に属する対称操作のループ
       do iopr_ind = 1, ks%kgroup%cls(ic)%order
!対称操作の名前の出力
          write(nfout, 403) '#', ks%kgroup%op(ks%kgroup%cls(ic)%ind(iopr_ind))%name
!可約表現の次元のループ
          do id = 1, dim
!同上
             do jd = 1, dim
!表現行列のid,jd成分の出力
                write(nfout, 402) '#', id, jd,&
                     redrep(ideg,is)%wd(id, jd, ks%kgroup%cls(ic)%ind(iopr_ind))
             end do
          end do
       end do
    end do
401 format(2(A,1X),'(',f0.8,',',f0.8,')')
402 format(A,1X,2(I0,1X),'(',f0.8,',',f0.8,')')
403 format(2(A,1X))
  end subroutine print_cr_wd
end module m_reduce_band
