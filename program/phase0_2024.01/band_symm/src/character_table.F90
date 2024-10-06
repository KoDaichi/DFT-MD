!指標表の処理
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_character_table
#include "defines.F"
  use m_commons, only: nearlyEquals, intcomp, complexguess
  use m_tspace_defines, only: irrep, character_table, group, kvector, num_imr_max, num_operator_max
  use m_tspace_commons, only: getKpointNameForTex
  use m_container_psicoef, only: ndim_spinor

  implicit none
  DOUBLE_T, parameter, private:: prec = 1e-5
#define NUM_MAX_CLASS_COLUMN 6
#define TEX_GROUP_NAME(SCHNAME) (SCHNAME(1:1)//'_{'//SCHNAME(2:3)//'}')
#define TEX_CHARTITLE(K,SCHNAME,HMNAME) (trim(K)//'('//TEX_GROUP_NAME(SCHNAME)//','//trim(HMNAME)//')')

contains

!行列のトレースをもとめる(表現行列の指標を求めるときに用いる。)
  COMPLEX_T function get_trace(wd, dim) result(cr)
    integer, intent(in):: dim!行列の次元
    COMPLEX_T, intent(in):: wd(6,6)!複素数型行列
    integer:: i
    cr = 0d0
    do i = 1, dim
       cr = cr + wd(i,i)
    end do
  end function get_trace

#include "get_alknown_ctable.F90"
#include "get_dgroup_ctable.F90"

!2つの表現a,bが同値である、つまり表現の指標が一致するとき.true.、同値でないとき.false.を返す
  logical function is_same_imr(a, b, order)
    type(irrep), intent(in):: a, b!表現構造体
    integer, intent(in):: order!2つの表現の類の数
    integer:: m
    COMPLEX_T:: c
    is_same_imr = .true.
    do m = 1, order
       c = a%comp(m) - b%comp(m)
       is_same_imr = is_same_imr .and. &
            nearlyEquals((real(c)**2 + aimag(c)**2),0d0,prec)
    end do

  end function is_same_imr

!mullikenの点群から既約表現を取得し、それと一致する表現にはmullikenの既約表現の名前をつける
  subroutine repname_remake(ctab, kgp)
    type(character_table), intent(inout):: ctab!TSPACEからもとめた既約表現
    type(group), intent(in):: kgp!群
    type(character_table):: ctab_new!mullikenの表現
    integer:: iimr, jimr
    logical:: status
!-----------------Modified by T.A.Ariasoca--------------------------------------
    if(ndim_spinor == 2) then
        ctab_new = get_dgroup_ctable(kgp%schname, status)
    else
        ctab_new = get_alknown_ctable(kgp%schname, status)!群のシェーンフリース記号からマリケンの既約表現を取得する
    end if
!-------------------------------------------------------------------------------
    if(.not.status) return
    do iimr = 1, ctab%numimr
       do jimr = 1, ctab%numimr
          if(is_same_imr(ctab%imrs(iimr), ctab_new%imrs(jimr),&
               ctab%nclass)) then
             ctab%imrs(iimr)%name_mulliken = ctab_new%imrs(jimr)%name!mullikenの既約表現と一致するとき、mullikenの既約表現の名前を保存する。
             cycle
          end if
       end do
    end do
  end subroutine repname_remake

!TSPACE内で記憶されているk点の群における既約表現の指標を返す。TSPACE関数tsirnr, irmataを使用
  type(character_table) function get_character_table(kp, kgroup)&
       result(ch)
#define NDIM 6
    type(kvector), intent(in):: kp!k点の名前、成分など
    type(group), intent(in):: kgroup!k点の群
    integer:: dime(num_imr_max), numimr, order, ic

    ch%maxdim = 1
!TSPACE内での指標表の設定
    call tsirnr(numimr, order, dime)
    ch%nclass = kgroup%nclass
    ch%numimr = numimr
    ch%numop = order
    call setchar(ch, 0)

!群構造体から指標表構造体へ、類の名前のコピー
    do ic = 1, kgroup%nclass
       ch%class_name(ic) = kgroup%cls(ic)%name
    end do
    call repname_remake(ch, kgroup)

  contains
    subroutine setchar(ch, istimr)
      type(character_table), intent(inout):: ch
      integer, intent(in):: istimr

      integer:: ic, i, j
      integer:: ir, mg, nd, jg(num_operator_max), jga(2, 3, num_operator_max) ! dummy
      COMPLEX_T:: cr(num_operator_max), wd(num_operator_max)
      DOUBLE_T:: dcr_inv

#define NAME(ir) ch%imrs(ir)%name
#define COMP(ir, ic) ch%imrs(ir)%comp(ic)
#define WD(i, j, ir, iopr) ch%imrs(ir)%wd(i, j, iopr)

    do ir = 1, numimr
!TSPACE内部データの読み取り
       call tsirmr(ir, 1, 1, mg, nd, jg, jga, cr, wd)
!指標表の構造体へコピー
       ch%imrs(ir + istimr)%dim = nd
       if(nd > ch%maxdim) ch%maxdim = nd

       do i = 1, nd
          do j = 1, nd
             call tsirme(i, j, wd)
             do ic = 1, ch%numop
                WD(i, j, ir + istimr, ic) = wd(ic)
             end do
          end do
       end do
       do ic = 1, ch%nclass
          COMP(ir + istimr,ic) = cr(kgroup%cls(ic)%ind(1))
       end do

!反転操作がある場合、その指標の二乗値を取得
       if(kgroup%inversion > 0) then
          dcr_inv = real(cr(kgroup%inversion))**2 + aimag(cr(kgroup%inversion))**2
       else
          dcr_inv = 0d0
       end if
!反転操作の指標から指標の名前を決定する
       if(.not.nearlyEquals(dcr_inv, 0d0, prec)) then
#define NUMB(ir) ((ir + 1)/2)
          if(mod(ir + istimr,2) == 0) then
             write(NAME(ir + istimr), '(A,I0,A)')&
                  '{',NUMB(ir + istimr), '}^-'
          else
             write(NAME(ir + istimr), '(A,I0,A)')&
                  '{',NUMB(ir + istimr), '}^+'
          end if
       else
          write(NAME(ir + istimr), '(A,I0,A)')&
               '{', ir + istimr, '}'
#undef NUMB
       end if
       ch%imrs(ir + istimr)%name_mulliken = ""
    enddo
  end subroutine setchar
end function get_character_table

!指標表の出力、既約表現と類にわけて出力される
  subroutine print_character_table(chtable, kgroup, nfout)
    type(character_table), intent(in):: chtable!指標表
    type(character_table):: ctab_new
    type(group), intent(in):: kgroup!k群
    integer, intent(in):: nfout!出力する装置番号
    integer:: ic, iimr, i, j
    logical:: status
    write(nfout, '(A)') '##Character'
    ctab_new = get_dgroup_ctable(kgroup%schname, status)
#define IMR chtable%imrs(iimr)
    do iimr = 1, chtable%numimr
       do ic = 1, chtable%nclass
          write(nfout, fmt = '(2A,1X,A24,1X)', advance = 'no')&
               '#', trim(IMR%name), trim(chtable%class_name(ic))
          write(nfout, 400) IMR%comp(ic)
          write(nfout, 400) ctab_new%imrs(iimr)%comp(ic)
       end do
       write(nfout, *)
    end do
    write(nfout, '(A)') '##Matrix'
    do iimr = 1, chtable%numimr
!       if(IMR%dim < 2) cycle
       do ic = 1, kgroup%order
          write(nfout, *) '#', trim(IMR%name), trim(kgroup%op(ic)%name), get_trace(IMR%wd(:,:,ic), IMR%dim)
          do i = 1, IMR%dim
             do j = 1, IMR%dim
                write(nfout, fmt = '(A,1X,2(I0,1X))', advance = 'no')&
                     '#', i, j
                write(nfout, 400) IMR%wd(i, j, ic)
             end do
          end do
          write(nfout, *)
       end do
    end do
#undef IMR
400 format('(',f0.8,',',f0.8,')')
  end subroutine print_character_table

!指標表をTeX形式で出力する
  subroutine print_character_table_for_TeX(cht, kgp, kp, nfout)
    type(character_table), intent(in):: cht!指標表
    type(group), intent(in):: kgp!k群
    type(kvector), intent(in):: kp!k点
    integer, intent(in):: nfout!出力する装置番号
#define WRITEA write(nfout,fmt='(a)',advance='no')

    write(nfout, *) TEX_BEGINTABLE
    write(nfout, *) TEX_ADJUST
!    write(nfout, *) TEX_CAPTION_FOR_CHARTABLE(trim(getKpointNameForTex(kp%name)))
    if(cht%nclass > NUM_MAX_CLASS_COLUMN) then
       write(nfout, '(A,I0,A)') TEX_BEGINTABULAR(cht%nclass/2+1)
       call print_table('$'//TEX_CHARTITLE(getKpointNameForTex(kp%name),kgp%schname,kgp%hmname)//'$',&
            1, cht%nclass/2)
       call print_table('', cht%nclass/2 + 1, cht%nclass)
    else
       write(nfout, '(A,I0,A)') TEX_BEGINTABULAR(cht%nclass+1)
       call print_table('$'//TEX_CHARTITLE(getKpointNameForTex(kp%name),kgp%schname,kgp%hmname)//'$',&
            1, cht%nclass)
    end if
    write(nfout, '(a)') TEX_ENDTABULAR
    call print_operator_in_class(nfout, kgp)
    write(nfout, *) TEX_ENDTABLE
    write(nfout, *)
  contains
!類と、類に属する対称操作の表の出力
    subroutine print_operator_in_class(nfout, kgp)
      integer, intent(in):: nfout
      type(group), intent(in):: kgp
      integer:: ic, iopr
      write(nfout, *) TEX_BEGINTABULAR(2)
      write(nfout, *) 'Class&Operator\\\hline'
      do ic = 1, kgp%nclass
         write(nfout, fmt = '(a)', advance = 'no')&
              trim(kgp%cls(ic)%name)//'&'//trim(kgp%op(kgp%cls(ic)%ind(1))%name)
         do iopr = 2, kgp%cls(ic)%order
            write(nfout, fmt = '(a)', advance = 'no')&
                 ','//trim(kgp%op(kgp%cls(ic)%ind(iopr))%name)
         end do
         write(nfout, *) '\\\hline'
      end do
      write(nfout, *) TEX_ENDTABULAR
    end subroutine print_operator_in_class
!類の一覧の出力
    subroutine print_class_table(isc, iec)
      integer, intent(in):: isc, iec
      integer:: ic
      do ic = isc, iec
         WRITEA '&$'//trim(cht%class_name(ic))//'$'
      end do
      write(nfout, '(a)') '\\ \hline'
    end subroutine print_class_table

!1つの既約表現の出力
!指標(複素数)を整数型複素数に変換し、その値から表示するデータを決める
    subroutine print_chars(imr, isc, iec)
      integer, intent(in):: imr, isc, iec
      type(intcomp):: icmp
      integer:: ic
      do ic = isc, iec
          WRITEA '&$'
          icmp = complexguess(cht%imrs(imr)%comp(ic))!整数型複素数の取得
          if(icmp%abs_type == 0) then!複素数の形が不明
10           write(nfout, fmt = 400, advance = 'no') cht%imrs(imr)%comp(ic)
             WRITEA '$'
             cycle
          end if

          if(icmp%abs_num == 0) then!複素数の形が不明
             WRITEA '0$'
             cycle
          end if

          if(icmp%arg_num == 0) then
             icmp%arg_type = 1 !arg = 0、つまり指標が正の実数のとき
          else if(icmp%arg_num == 1 .and. icmp%arg_den == 2 &
               .and. icmp%arg_type == 1) then
             icmp%arg_type = 2 !arg = pi/2、つまり指標が純虚数のとき
          else if(icmp%arg_num == 1 .and. icmp%arg_den == 2 &
               .and. icmp%arg_type == 2) then
             icmp%arg_type = 4 !arg = -pi/2、つまり指数が純虚数で、負の符号がつくとき
          else if(icmp%arg_num == 1 .and. icmp%arg_den == 1 &
               ) then
             icmp%arg_type = 3 !arg = pi、つまり指標が負の実数のとき
          else
             icmp%arg_type = icmp%arg_type*10
          end if

          if(icmp%abs_den == 1) then !ノルムが分数でない、つまり|z|=m or sqrt(m)
             icmp%abs_type = icmp%abs_type * 10 + 1
          else
             icmp%abs_type = icmp%abs_type * 10 + 2
          end if

!指標のノルムが整数のとき(abs_type=11)、平方根のとき(abs_type==21)、分数のとき(abs_type==12)、平方根/整数のとき(abs_type=22)に合わせて指標のノルムだけを出力するマクロ、ただしノルムが1のときは出力する必要がないので、そのときは出力しない
#define ABSGUESS\
          select case(icmp%abs_type);\
          case (11);\
             if(icmp%abs_num /= 1) write(nfout, fmt = '(I0)', advance = 'no') icmp%abs_num;\
          case (21);\
             write(nfout, fmt = '(A,I0,A)', advance = 'no') '\sqrt{', icmp%abs_num, '}';\
          case (12);\
             write(nfout, fmt = '(A,I0,A,I0,A)', advance = 'no') '\frac{', icmp%abs_num, '}{', icmp%abs_den, '}';\
          case (22);\
             write(nfout, fmt = '(A,I0,A,I0,A)', advance = 'no') '\frac{\sqrt{', icmp%abs_num, '}}{', icmp%abs_den, '}';\
          end select
!マクロ定義終わり


          select case(icmp%arg_type)
          case(1)!arg = 0、指標が正の実数のとき
             ABSGUESS
             if(icmp%abs_num == 1) WRITEA '1'
          case(2)!arg = pi/2、指標が正の純虚数のとき
             ABSGUESS
             WRITEA 'i'
          case(3)!arg = pi、指標が負の実数のとき
             WRITEA '-'
             ABSGUESS
             if(icmp%abs_num == 1) WRITEA '1'
          case(4)!arg = -pi/2、指標が負の純虚数のとき
             WRITEA '-'
             ABSGUESS
             WRITEA 'i'
          case(10)!指標が複素数（純虚数でない）かつ位相が正のとき
             ABSGUESS
             write(nfout, fmt = '(A,I0,A,I0,A)', advance = 'no') 'e^{i\frac{', icmp%arg_num, '}{', icmp%arg_den, '}\pi}'
          case(20)!指標が複素数（純虚数でない）かつ位相が負のとき
             ABSGUESS
             write(nfout, fmt = '(A,I0,A,I0,A)', advance = 'no') 'e^{-i\frac{', icmp%arg_num, '}{', icmp%arg_den, '}\pi}'
#undef ABSGUESS
          case default
             go to 10
          end select
          WRITEA '$'
       end do
400 format('(',f0.8,',',f0.8,')')
     end subroutine print_chars

!指標の出力
     subroutine print_table(title, isc, iec)
       character(*), intent(in):: title
       integer, intent(in):: isc, iec
       integer:: imr

       write(nfout, fmt = '(A)', advance = 'no') title
       call print_class_table(isc, iec)

       do imr = 1, cht%numimr
          write(nfout, fmt = '(A)', advance = 'no') '$'//trim(cht%imrs(imr)%name)//'$'
          if(cht%imrs(imr)%name_mulliken /= "")&
               write(nfout, fmt = '(A)', advance = 'no') '($'//trim(cht%imrs(imr)%name_mulliken)//'$)'
          call print_chars(imr, isc, iec)
          write(nfout, '(a)') '\\ \hline'
       end do
     end subroutine print_table
  end subroutine print_character_table_for_TeX
end module m_character_table
