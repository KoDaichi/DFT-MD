!k点の適合関係についての処理
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_compatibility
#include "defines.F"
  use m_commons, only: stdout
  use m_tspace_defines, only: kvector, character_table, comtable
  use m_tspace_commons, only: getKpointNameForTex
!  use m_character_table
  implicit none

!k点バッファの最大要素数
  integer, parameter, private:: num_kpoint_max = 255
!k点バッファの現在の要素数
  integer, private:: numk = 0
!k点のバッファ(k-pathのすべての点ではなく、k-pathの各節と中間から1つずつ登録する
  type(kvector), private:: kps(num_kpoint_max)
  type(character_table), private:: chts(num_kpoint_max)
contains 
!k点を二つ渡し(ik1, ik2)、それら二つのk点群間の適合関係を返す。ik1がik2の部分群になるように指定する。isdoubleは任意変数、TRUEにすると二重表現における適合関係を返す。errは計算が成功したかどうか、部分群の条件を満たさないなどのエラーが発生したときFALSEが格納される。TSPACE関数compatを使用
  type(comtable) function get_compat_table(ik1, ik2, err) result(tbl)
    type(kvector), intent(in):: ik1, ik2
    logical, intent(out):: err
    integer:: isd, ind, jsd
    integer:: ndes(12), jtrs(12), ipas(12)
    tbl%kp(1) = ik1
    tbl%kp(2) = ik2
    if(ik1%isdouble) then
       isd = 1
    else
       isd = 0
    endif

    tbl%comp = 0
    call compat(ik1%comp, ik1%common_denom, tbl%nimr(1),&
         ik2%comp, ik2%common_denom, tbl%nimr(2),&
         isd, tbl%comp, ind)
    !call cmptrv(ik1%comp, ik1%common_denom, ik2%comp, ik2%common_denom, isd, tbl%comp, tbl%comp)
    if(ind > 0) then
       err = .false.
    else
       err = .true.
    endif
    !call dgtrst(jdub,nrs,mmg,nstr,ndes,jtrs,ipas)

  end function get_compat_table

!k点のリストkpsに入力されたkpが登録されていたら.true.、されていなかったら.false.を返す
!k点の比較はk点の名前同士を比較する
  logical function isDuplicateKpoint(kp) result(is)
    type(kvector), intent(in):: kp  
    integer:: ik
    is = .false.
    do ik = 1, numk
       if(kps(ik)%name(1:2) == kp%name(1:2)) is = .true.
    end do
  end function isDuplicateKpoint

!k点バッファにkpを追加する。重複チェックは行われない。バッファの数が最大数を超えた時、バッファ先頭の値が上書きされる
 subroutine add_kpoint(kp, cht)
   type(kvector), intent(in):: kp
   type(character_table), intent(in):: cht
!k点の追加
   numk = numk + 1
   numk = mod(numk, num_kpoint_max)
   kps(numk) = kp
   chts(numk) = cht
 end subroutine add_kpoint

!k点のリストkpsを、k点のインデックスをつかってソートする
 subroutine sort_kpoint()
   integer:: i, j
   type(kvector):: kp
   type(character_table):: cht
   j = numk - 1
   do while(j > 1)
      do i = 1, j
         if(kps(i)%index > kps(i + 1)%index) then
            kp = kps(i)
            kps(i) = kps(i + 1)
            kps(i + 1) = kp

            cht = chts(i)
            chts(i) = chts(i + 1)
            chts(i + 1) = cht
         end if
      end do
      j = j - 1
   end do
 end subroutine sort_kpoint

!適合関係出力のためのフォーマット
!あるk点(knameとは無関係)のある既約表現imr_sがk点(knameと対応する)のどの既約表現と適合するかを、適合関係comtを用いて判断し、出力を行う。
  subroutine format_compat_table(imr_s, kname, cht, comt, nfout)
    integer, intent(in):: imr_s, nfout 
    character(len=8), intent(in):: kname
    type(character_table), intent(in):: cht
    type(comtable), intent(in):: comt
    integer:: imr_l, imr_k!k点(knameと対応する)の既約表現のインデックス
#define COMP(i) comt%comp(imr_s, i)
    imr_k = 1
!TSPACE関数compatから得た、既約表現imr_sの中に部分群の既約表現imr_lが含まれる数をcomt%comp(imr_s, imr_l)が記憶している
!最初に含まれる既約表現の番号までスキップし、その番号をimr_kに記憶させる
    do while(COMP(imr_k) == 0)
       imr_k = imr_k + 1
       if(imr_k > comt%nimr(1)) go to 200
    end do
    if(COMP(imr_k) > 1) write(nfout, fmt = '(I0)', advance = 'no') COMP(imr_k)
    write(nfout, fmt = '(A)', advance = 'no') TEX_COMPAT_NAME(kname, cht%imrs(imr_k)%name)
!最初以降の既約表現が含まれるかどうかを調べる
    do imr_l = imr_k + 1, comt%nimr(1)
       select case(COMP(imr_l))
       case(0)!1つも含まれない場合
          cycle
       case(1:)!1つ以上含まれる場合
          if(COMP(imr_l) > 1) write(nfout, fmt = '(I0)', advance = 'no') COMP(imr_l)
          write(nfout, fmt = '(2A)', advance = 'no') '+',TEX_COMPAT_NAME(kname, cht%imrs(imr_l)%name)
       end select
    end do
200  continue
#undef COMP
  end subroutine format_compat_table

!適合関係の出力(バッファ内にあるk点から2つ選ぶ、すべてのk点の組み合わせについて適合関係を出力する)
  subroutine print_compat_table(nfout)
    integer, intent(in):: nfout
    integer:: ik

    call sort_kpoint
!    open(stdout)
    do ik = 1, numk - 1
       call PRINTCT(nfout, ik, ik + 1)
       call PRINTCT(nfout, ik + 1, ik)
    end do
!    close(stdout)
    return
contains
!適合関係の出力関数(k点を2つ指定する)
  subroutine PRINTCT(nfout, ik1, ik2)
    integer, intent(in):: nfout, ik1, ik2
    type(comtable):: comt
    integer:: imr_s
    logical:: is_success
!適合関係の取得
    comt = get_compat_table(kps(ik1), kps(ik2), is_success)
    if(.not.is_success) then
       write(stdout, *) '!COMTABLEERROR', ik1, kps(ik1)%name, ik2, kps(ik2)%name
       return
    end if
!TeXのtabular環境の出力
    write(nfout, *) TEX_BEGINTABULAR(2)
    write(nfout, *) '$', getKpointNameForTex(kps(ik2)%name), '$&$',&
         getKpointNameForTex(kps(ik1)%name), '$ \\ \hline\hline'
    do imr_s = 1, comt%nimr(2)
       write(nfout, fmt = '(3A)', advance = 'no') '$',TEX_COMPAT_NAME(getKpointNameForTex(kps(ik2)%name),chts(ik2)%imrs(imr_s)%name),'$&$'
       call format_compat_table(imr_s, getKpointNameForTex(kps(ik1)%name), chts(ik1), comt, nfout)
       write(nfout, '(A)') '$\\ \hline'
    end do
    write(nfout, *) TEX_ENDTABULAR
  end subroutine PRINTCT
  end subroutine print_compat_table

  subroutine print_compat_csv(nfout)
    integer, intent(in):: nfout
    integer:: ik

!    open(stdout)
    call sort_kpoint
    do ik = 1, numk - 1
       call PRINTCT(nfout, ik, ik + 1)
       call PRINTCT(nfout, ik + 1, ik)
    end do
!    close(stdout)
    return
contains
!適合関係の出力関数(k点を2つ指定する)
  subroutine PRINTCT(nfout, ik1, ik2)
    integer, intent(in):: nfout, ik1, ik2
    type(comtable):: comt
    integer:: imr_s
    logical:: is_success
    comt = get_compat_table(kps(ik1), kps(ik2), is_success)
    if(.not.is_success) then
       write(stdout, *) '!COMTABLEERROR', ik1, kps(ik1)%name, ik2, kps(ik2)%name
       return
    end if
    write(nfout, *) getKpointNameForTex(kps(ik2)%name)//','//getKpointNameForTex(kps(ik1)%name)
    do imr_s = 1, comt%nimr(2)
       write(nfout, fmt = '(A)', advance = 'no') TEX_COMPAT_NAME(getKpointNameForTex(kps(ik2)%name),chts(ik2)%imrs(imr_s)%name)//','
       call format_compat_table(imr_s, getKpointNameForTex(kps(ik1)%name), chts(ik1), comt, nfout)
       write(nfout, *)
    end do
    write(nfout, *)
  end subroutine PRINTCT

  end subroutine print_compat_csv

end module m_compatibility
