!TSPACEを使うためのユーティリティ、以下は変数の定義部分
!Copyright (C) 2014-2015 Ryosuke Tomita 

#define DEBUG(F) write(100,F)
module m_tspace_commons
#include "defines.F"
!  use m_character_table, only: get_trace
  use m_commons, only: pi, stdout
  use m_phase_commons, only: get_point_group_name
  use m_tspace_defines, only: num_operator_max, class, kvector, group, lattice_parameter, num_class_max, num_operator_per_class_max, element_name_inversion, element_name_cubic, element_name_hexagonal
  use m_container_psicoef,only:ndim_spinor ! Saito
  implicit none

#define IOPR_TSP(i) ig(jg(i))
!TSPACEの内部変数
!(以下、OhもしくはD6h群での番号)対称操作のJones表現、群表(O_aO_b == O_im(a,b))、逆元表
  integer, private:: it(nxyz_, num_operator_max), im(num_operator_max, num_operator_max), iv(num_operator_max)
!ブラベ格子の番号、空間群の位数、(以下、OhもしくはD6h群での番号)空間群の対称操作の番号、空間群の並進ベクトル
  integer, private:: il, ng, ig(num_operator_max), jv(2, nxyz_, num_operator_max)
!k点の分母、共通分子、k点群の位数、(以下、Γ点群での番号)k点群の対称操作の番号、回転操作に付随するk点に等価な点の番号
  integer, private:: kb(nxyz_), icb, mg, jg(num_operator_max), jk(nxyz_, num_operator_max)
!スピン空間の回転行列(スピン状態は反転対称なので位数の最大は24)
  COMPLEX_T, private:: sn(2,2,num_operator_max/2)
  common/spg1/it, im, iv
  common/spg2/il, ng, ig, jv
  common/spg3/kb, icb, mg, jg, jk
  common/spg6/sn

contains
#ifdef DEBUG
  subroutine debugjv(nfout)
    integer, intent(in):: nfout
    integer:: iopr, j
    write(nfout, *) ig
    do iopr = 1, ng
       write(nfout, *) iopr, ig(iopr), (jv(1:2, j, iopr), j = 1, 3)
    end do
  end subroutine debugjv
#endif

  subroutine tescom()
!      type(character_table), pointer, private:: chtable
      integer:: ndes(12), jtrs(12),ipas(12)
      integer:: i, jdub, nrs, mmg, nstr
!        i= chtable%numimr
      call dgtrst(jdub, nrs, mmg, nstr, ndes, jtrs, ipas)
!      write(stdout,*)'order of kpoint group =', mmg, ' number of star =', nstr,' number of representation =', nrs
      do i = 1, nrs
          if(jtrs(i)==0) write(stdout,*)  '     PAIRING of ',i,' and',ipas(i),':    Herring sum=',jtrs(i)
          if(jtrs(i) < 0) write(stdout,*) '     DUBLING of ',i, ':    Herring sum=',jtrs(i)
!          write(stdout, *)'no', i
!          write(stdout, *)'degeneracy = ', ndes(i)
!          write(stdout, *)'herring sum = ', jtrs(i)
!          write(stdout, *)'partner no = ', ipas(i)
      end do

  end subroutine tescom
!対称操作のインデックス及びTSPACE内部変数から対称操作の文字列表現を得る
  character(len = NameLen) function get_element_name(iopr)&
       result(name)
    integer, intent(in):: iopr
    name = ''
    if(il /= 0) then !立体群(Oh)
       if(iopr > NOP_MAX_CUB) then
          name = ''
!番号が最大値の半分を超えるとき、文字列先頭に'I'をつける
       else if(iopr > NOP_MAX_CUB/2) then
          name = trim(element_name_inversion) // trim(element_name_cubic(iopr - NOP_MAX_CUB/2))
       else
          name = trim(element_name_cubic(iopr))
       end if
    else !六方晶群(D6h)
       if(iopr > NOP_MAX_HEX) then
          name = ''
       else if(iopr > NOP_MAX_HEX/2) then
          name = trim(element_name_inversion) // trim(element_name_hexagonal(iopr - NOP_MAX_HEX/2))
       else
          name = trim(element_name_hexagonal(iopr))
       end if
    end if
  end function get_element_name

!類に属する対称操作の名前から、類の名前を決定する(beta版)
  character(len = NameLen) function get_class_name(cls,grp)&
       result(name)
    type(class), intent(in):: cls
    type(group), intent(in):: grp
    if(cls%order == 1) then
       name = trim(grp%op(cls%ind(1))%name)
    else
       write(name, '(I0,A)') cls%order, trim(grp%op(cls%ind(1))%name)
    end if
  end function get_class_name

!k点群の対称操作の共役元をもとめ、類に分ける,出力は類の数
  subroutine set_class(kpoint, kgroup)
#define CLASS kgroup%cls
#define ORDER kgroup%order
#define NCLASS kgroup%nclass
    type(kvector), intent(in):: kpoint
    type(group), intent(inout):: kgroup
    integer:: iop1, iop2, iop_, jg1, jg2, imbuf
!ある対称操作が類に分けられたかどうか
    logical:: iscount(num_operator_max)
    integer:: v1(3), v2(3), maxang

    NCLASS = 0
    maxang = 0
    
    iscount = .false.
    do iop1 = 1, ORDER!iterate through number of order of kgroup  
       if(NCLASS >= num_class_max) then
          NCLASS = num_class_max
          exit
       end if
!対称操作がすでに処理されているとき、スキップ
       if(iscount(iop1)) cycle
       NCLASS = NCLASS + 1!this iteration set the number of class
       CLASS(NCLASS)%order = 1
       CLASS(NCLASS)%ind = 0
       CLASS(NCLASS)%ind(CLASS(NCLASS)%order) = iop1
       iscount(iop1) = .true.
!対称操作のインデックス(TSPACE内でのインデックス)
       jg1 = IOPR_TSP(iop1)!A
       do iop2 = 1, ORDER
          if(iop2 == iop1) cycle
          jg2 = IOPR_TSP(iop2)!B
!write(stdout, *) iop1, jg1, iop2, jg2, NCLASS, CLASS(NCLASS)%order
!群表imおよび逆元ivにより共役元BAB^-1のインデックスをもとめる
          imbuf = im(im(jg2, jg1), iv(jg2))!BAB^-1
          do iop_ = 1, ORDER!C
!共役元がk点群の中にあり、まだ処理されていないとき
             if(.not.iscount(iop_) .and. imbuf == IOPR_TSP(iop_))& 
                 go to 10!C == BAB^-1           
          end do
          cycle
!類表に共役元を追加する
10        if(CLASS(NCLASS)%order <= num_operator_per_class_max)&
               CLASS(NCLASS)%order = CLASS(NCLASS)%order + 1
          !write(stdout,*)iop1,iop2,iop_,'|',CLASS%order 
          CLASS(NCLASS)%ind(CLASS(NCLASS)%order) = iop_
!この対称操作は処理したのでフラグを立てる
       iscount(iop_) = .true.
       end do
    end do
!類の名前付け
    do jg1 = 1, NCLASS
       CLASS(jg1)%name = get_class_name(CLASS(jg1), kgroup)
    end do

#undef CLASS
#undef ORDER
#undef NCLASS
  end subroutine set_class

!対称操作から点群の名前を得る
  subroutine set_pg_name(grp)
    type(group), intent(out):: grp
    integer:: i, nopr(num_operator_max)
    grp%schname = ''
    grp%hmname = ''
    do i = 1, mg
       nopr(i) = IOPR_TSP(i)
!write(stdout, *) i, nopr(i)
    enddo
    if(il /= 0) then
       call get_point_group_name(stdout, 'cubic    ',&
            mg, nopr, grp%schname, grp%hmname)
    else
       call get_point_group_name(stdout, 'hexagonal',&
            mg, nopr, grp%schname, grp%hmname)
    end if
  end subroutine set_pg_name

!Γ点群を渡し、TSPACEで指定されているk点の群を得る
  type(group) function get_kgroup_from_gammagroup(kgp)&
       result(buf)
    type(group), intent(in):: kgp
    integer:: i
    buf%order = mg
    do i = 1, buf%order
       buf%op(i) = kgp%op(jg(i))
    enddo
  end function get_kgroup_from_gammagroup

!TSPACEから空間群もしくはk点群を得る関数
  type(group) function get_group_template(is_kgroup) result(gr)
    logical, intent(in):: is_kgroup
    integer:: i, ittmp, j, k, igshift, max_
!k点群かどうかで位数が変わる
!六方晶のk点群の場合, 回転軸が実格子と異なるための措置、詳しくはTSPACEを参照されたし
    igshift = 0
    if(is_kgroup) then
       gr%order = mg!this set kgroup order, mg taken from TSPACE
       if(il == 0)&
            igshift = num_operator_max/2
    else
       gr%order = ng
    end if
    if(il == 0) then
       max_ = num_operator_max/4
    else
       max_ = num_operator_max/2
    end if

    gr%inversion = 0
    do i = 1, gr%order
       gr%op(i)%index_for_gamma = jg(i)
       gr%op(i)%index_for_all = IOPR_TSP(i)

       if((il == 0 .and. IOPR_TSP(i) == 13) .or.&
            (il /= 0 .and. IOPR_TSP(i) == 25))&
            gr%inversion = i
       gr%op(i)%rot = 0d0
       gr%op(i)%tr = 0d0
!k点群は空間群の回転操作のみからなる部分群である。k点群においてどの回転操作が要素となるかはk点によって異なり、それについての情報がTSPACE内部変数jg(1: mg)に記憶されている
       do j = 1, nxyz_
!TSPACE内部変数itに対称操作の回転についての情報が格納されている。jones表現XYZはそれぞれ整数1,2,3で表される。-X, -Y, -Z は-1, -2, -3として表現される
          ittmp = it(j, igshift + IOPR_TSP(i))

!jones表現X-Y(TSPACE内でWと表現される)は整数4として記憶されている。六方晶k点群の場合、Wの定義はX+Yに変わる。符号の反転はXYZと同様
          if(abs(ittmp) == 4) then
             k = sign(1, ittmp)
             gr%op(i)%rot(1, j) = k
             gr%op(i)%rot(2, j) = -k
             if(is_kgroup) gr%op(i)%rot(2, j) = k
          else
             k = abs(ittmp)
             gr%op(i)%rot(k, j) = sign(1, ittmp)
          endif
!並進についての情報がTSPACE内部変数jvに記憶されている。jv(1, 1:3, 1:48)が並進ベクトルの分子、jv(2, 1:3, 1:48)が分母、ifは分母が0であったときの措置
          if(jv(2, j, jg(i)) == 0) then
             gr%op(i)%tr(j) = 0d0
          else
             gr%op(i)%tr(j) = dble(jv(1, j, jg(i)))&
                  /jv(2, j, jg(i))
          end if
       enddo
#ifdef _USE_LAPACK_
       gr%op(i)%vaxis = get_rotation_axis(gr%op(i)%rot)
       gr%op(i)%vang = .5d0*(gr%op(i)%rot(1,1) + gr%op(i)%rot(2,2) + gr%op(i)%rot(3,3) - 1d0)
#endif
       gr%op(i)%spinrot(:,:) = sn(:,:,mod(IOPR_TSP(i) - 1, max_) + 1)
       gr%op(i)%name = trim(get_element_name(IOPR_TSP(i)))
    enddo
  end function get_group_template

!TSPACE内で(TSPACE関数TSIREPにより)設定されているk点についての対称操作の情報を返す
  type(group) function get_kgroup() result(gr)
    gr = get_group_template(.true.)
  end function get_kgroup

!TSPACEで設定された結晶系(TSPACE関数TSPACEなどで設定)の空間群の対称操作の情報を返す
  type(group) function get_group() result(gr)
    gr = get_group_template(.false.)
  end function get_group

!晶系を表す文字列[FICRHP]から、TSPACEに指定する整数ilを決定する(inputを分かり易くするための措置)
  integer function determine_lattice_system(hmname) result(il)
    character(len = NameLen), intent(in):: hmname
    select case(hmname(1:1))
    case ('f', 'F')
       il = 2
    case ('i', 'I')
       il = 3
    case ('c', 'C')
       il = 4
    case ('r', 'R')
       il = -1
    case ('h', 'H')
       il = 0
    case ('p', 'P')
       il = 1
    case default
       il = 255
    end select
  end function determine_lattice_system

!結晶系を決める整数(上記のdetermine_lattice_systemから求められる)と生成元および反転対称性の有無から、TSPACE内で空間群を決定するための関数(TSPACE関数tspace, tsgenr, tspgrpを使用)
  subroutine determine_space_group(latsys, ngen, noprs, ntrs, inv)
    integer, intent(in):: latsys, ngen, noprs(ngen), ntrs(2, nxyz_, ngen)
    logical, intent(in):: inv
    integer:: i
!結晶系の決定
    call tspace(latsys)
    do i = 1, ngen
!生成元の登録
       call tsgenr(noprs(i), ntrs(:, :, i))
    enddo
!反転対称性の有無
    if(inv) then
       call tspgrp(1)
    else
       call tspgrp(0)
    endif
  end subroutine determine_space_group

!TSPACE内で格子定数を決定する(TSPACE関数tslatcを使用)
  subroutine set_lattice_parameters(lp)
    type(lattice_parameter), intent(in):: lp
!結晶系の設定
    call determine_space_group(determine_lattice_system(lp%lattice_type),&
            lp%num_generator, lp%index_rotation, &
            lp%translation, lp%isinversion)
!格子定数の設定
    call tslatc(lp%a, lp%b, lp%c,&
         cos(lp%alpha * pi/180),&
         cos(lp%beta * pi/180),&
         cos(lp%gamma * pi/180))
  end subroutine set_lattice_parameters

!TSPACE内でk点を指定する(TSPACE関数tsirep, kpnameを使用)
  subroutine set_kpoint_irep(kp, nspin)
    type(kvector), intent(inout):: kp
    integer, intent(in), optional:: nspin
    if(present(nspin) .and. nspin < 2) then
       call tsirep(kp%comp, kp%common_denom, 0)
!    else if(kp%isdouble) then
    else if(ndim_spinor == 2) then
       call tsirep(kp%comp, kp%common_denom, 1)
    else
       call tsirep(kp%comp, kp%common_denom, 0)
    end if
    kp%name = ""
    call kpname(kp%comp(1), kp%comp(2), kp%comp(3), kp%common_denom, kp%name)
    call tskfbz(kp%comp(1), kp%comp(2), kp%comp(3), kp%common_denom, kp%comm2bz)
  end subroutine set_kpoint_irep

  subroutine print_kpoint_and_kgroup(nfout, kp, grp)
    integer, intent(in):: nfout
    type(kvector), intent(in):: kp
    type(group), intent(in):: grp
    character*2 ctempo
    data ctempo/'  '/ 
    write(nfout, '(I0,1X,A,1X,4(I0,1X),2(A,1X),I0)') kp%index, kp%name, kp%comp, kp%common_denom, grp%schname, grp%hmname, kp%comm2bz
!     PRINT HERRING SUM  Saito
     if (ctempo.ne.kp%name)  call tescom()
     ctempo=kp%name 
!     
  end subroutine print_kpoint_and_kgroup    

!群の対称操作の出力
  subroutine print_group(nfout, grp)
    integer, intent(in):: nfout
    type(group), intent(in):: grp
    integer:: iopr, jd
#define OP grp%op(iopr)
    write(nfout, '("order=",i6)')grp%order
    do iopr = 1, grp%order
       write(nfout, *) iopr, IOPR_TSP(iopr), OP%name, OP%vang
       do jd = 1, 3
          write(nfout, *) OP%rot(:, jd), OP%tr(jd), OP%vaxis(jd)
       end do
       write(nfout, *) 'SpinRot'
       do jd = 1, 2
          write(nfout, *) OP%spinrot(jd,:)
       end do
    end do
#undef OP
  end subroutine print_group

!類の出力、フォーマット=> (類の位数、１のとき省略)(類の名前):(類に属する対称操作１の名前) (対称操作２の名前) ...
  subroutine print_class(nfout, kgroup)
    type(group), intent(in):: kgroup
    integer, intent(in):: nfout
    integer:: ic, iop
    write(nfout, '(a)') '#OPERATOR IN CLASS'
!write(nfout, *) kgroup%cn, kgroup%sv, kgroup%sh, kgroup%inv
    do ic = 1, kgroup%nclass
       write(nfout, fmt = '(A,1X,A10,2X)', advance = 'no')&
            '#class', trim(kgroup%cls(ic)%name)//':'
       do iop = 1, kgroup%cls(ic)%order
          write(nfout, fmt = '(A,1X)', advance = 'no')&
               trim(kgroup%op(kgroup%cls(ic)%ind(iop))%name)
       end do
       write(nfout, *)
    end do
  end subroutine print_class

  character(len=8) function getKpointNameForTex(name1) result(name2)
    character(len=2), intent(in):: name1
    select case(name1(1:2))
    case("GM")
       name2 = "\Gamma"
    case("LD")
       name2 = "\Lambda"
    case("SM")
       name2 = "\Sigma"
    case("DT")
       name2 = "\Delta"
    case("TP")
       name2 = "T'"
    case("SP")
       name2 = "P'"
    case default
       name2 = name1
    end select
  end function getKpointNameForTex
end module m_tspace_commons
