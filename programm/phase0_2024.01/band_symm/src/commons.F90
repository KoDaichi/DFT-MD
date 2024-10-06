!物理定数、数学公式など
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_commons
#include "defines.F"
  implicit none
!標準入力、標準出力、標準エラー出力の装置番号
  integer, parameter:: stdin = 5, stdout = 6, stderr = 0
!pi,2pi
  DOUBLE_T, parameter:: pi = acos(-1d0), pi2 = 2*pi
!3x3単位行列
  DOUBLE_T, dimension(3, 3), parameter:: &
       unit_matrix = reshape((/1., 0., 0., &
       0., 1., 0., &
       0., 0., 1./), (/3, 3/))
!複素単位
  COMPLEX_T, parameter:: zi = cmplx(0,1)
!複素数の比較(complexguess)に用いる、整数型複素数
  type intcomp
!複素数のノルム部分|z|の整数型。type==1のとき、|z|=num/den,type==2のとき、|z|=sqrt(num)/denとなるように設定する
     integer:: abs_num, abs_den, abs_type
!複素数の位相部分t=arctan(z/|z|)の整数型。type==1のとき、t=num/den*pi,type==2のとき、t=-num/den*piとなるように設定する
     integer:: arg_num, arg_den, arg_type
  end type intcomp

contains
!引き渡したファイルのパスを用いてファイルを読み込みモードで開き、引き渡した装置番号と関連付ける。開けた場合は.true.、開けなかった場合は.false.を返す。
  logical function openInput(fileno, path) result(is)
    integer, intent(in):: fileno!ファイルの装置番号
    character(*), intent(in):: path!ファイルのパス
    is = .false.
    open(fileno, file = trim(path), status = 'old', err = 100)
    is = .true.
100 return
  end function openInput

!コマンドライン引数をファイルのパスとみなし、そのファイルを開く。開けた場合は.true.、開けなかった場合は.false.を返す。
  logical function openInputFromArg(fileno) result(is)
    integer, intent(in):: fileno!ファイルの装置番号
    character(len=255):: argbuf
    integer, save:: nc = 1
    is = .false.
    if(nc > iargc()) return
    call getarg(nc, argbuf)
    nc = nc + 1
    if(.not.openInput(fileno, argbuf)) return
    is = .true.
  end function openInputFromArg

!コマンドラインからファイルを読み込みモードで開くが、失敗した時は引数のパスからファイルを読み込みモードで開く。開けた場合は.true.、開けなかった場合は.false.を返す。
  logical function openInputFromArgOrDefaultPath(fileno, path) result(is)
    integer, intent(in):: fileno!ファイルの装置番号
    character(*), intent(in):: path!コマンドラインからのファイル読み込みが失敗したときに読み込まれるファイル
    is = openInputFromArg(fileno) .or. openInput(fileno, path)
  end function openInputFromArgOrDefaultPath

!引き渡したファイルのパスを用いてファイルを書き込みモードで開き、引き渡した装置番号と関連付ける。開けた場合は.true.、開けなかった場合は.false.を返す。
  logical function openOutput(fileno, path)
    integer, intent(in):: fileno!ファイルの装置番号
    character(*), intent(in):: path!ファイルのパス
    openOutput = .false.
    open(fileno, file = trim(path), status = 'replace', err = 101)
    openOutput = .true.
101 return
  end function openOutput

!整数のベクトルが0である場合は.true.、0出ない場合は.false.を返す
  logical function is_vector_zero(kv)
    integer, dimension(3), intent(in):: kv!3次元整数ベクトル
    is_vector_zero = kv(1) == 0 .and. &
         kv(2) == 0 .and. &
         kv(3) == 0
  end function is_vector_zero

!2つの整数ベクトルが一致するとき.true.,一致しないとき.false.を返す
  logical function equals_vector(kv1, kv2)
    integer, dimension(3), intent(in):: kv1, kv2!3次元整数ベクトル1,2
    equals_vector = kv1(1) == kv2(1) .and. &
         kv1(2) == kv2(2) .and. &
         kv1(3) == kv2(3)
  end function equals_vector

!2つの3次元ベクトルの外積を計算して、3次元配列として返す
  function cross_product(a, b) result(c)
    DOUBLE_T, dimension(3), intent(in):: a, b!3次元ベクトルa,b
    DOUBLE_T, dimension(3):: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

!ベクトルのノルムを実数値として返す
  DOUBLE_T function norm(vect)
    DOUBLE_T, dimension(3), intent(in):: vect!ベクトル
    norm = sqrt(dot_product(vect, vect))
  end function norm

!引き渡されたベクトルを規格化する
  subroutine normalize(vect)
    DOUBLE_T, dimension(3), intent(inout):: vect!ベクトル
    vect = vect / norm(vect)
  end subroutine normalize

!格子パラメータα、β、γから規格化された格子ベクトルを求める。ただし、格子ベクトルv_1をx軸に並行であり、v_2はxy平面に存在するものとする。
!v_1 = e_x
!v_2 = cos(γ)e_x + sin(γ)e_y
!v_3 = cos(α)we_x + (cos(β)-cos(α)cos(γ))/sin(γ)e_y + g e_z
!g = sqrt(1 - cos(α)^2 - ((cos(β)-cos(α)cos(γ))/sin(γ))^2)
  function euler2vector(alpha, beta, gamma) result(v)
    DOUBLE_T, intent(in):: alpha, beta, gamma!格子パラメータ
    DOUBLE_T, dimension(3, 3):: v
    DOUBLE_T:: csc, snc, x, y, z
    csc = cos(gamma)
    snc = sin(gamma)
    v(:, 1) = (/1d0, 0d0, 0d0/)
    v(:, 2) = (/csc, snc, 0d0/)
    x = cos(alpha)
    y = (cos(beta) - x * csc)/snc
    z = sqrt(1 - x * x - y * y)
    v(:, 3) = (/x, y, z/)
  end function euler2vector

!実格子ベクトルa,b,cから逆格子ベクトルa~,b~,c~をもとめる。
!実格子ベクトルは3x3の行列として入力する。A(1:3,1)がa,A(1:3,2)がb,A(1:3,3)がcを意味する。
!逆格子ベクトルも3x3の行列として返す。。B(1:3,1)がa~,B(1:3,2)がb~,B(1:3,3)がc~を意味する。
  function get_recip_vector(A) result(B)
    DOUBLE_T, dimension(3, 3), intent(in):: A!実格子ベクトル
    DOUBLE_T, dimension(3, 3):: B!逆格子ベクトル
    DOUBLE_T:: c!実格子の体積
    B(:, 1) = cross_product(A(:, 2), A(:, 3))
    c = pi2/dot_product(A(:, 1), A(:, 1))
    B(:, 1) = b(:, 1) * c
    B(:, 2) = c * cross_product(A(:, 3), A(:, 1))
    B(:, 3) = c * cross_product(A(:, 1), A(:, 2))
  end function get_recip_vector

!実数値の比較。誤差errを含めてa,bが等しいときに.true.、異なるときに.false.を返す
!誤差errを含めてa,bが等しいとき|a-b|<errが成り立つ。
  logical function nearlyEquals(a,b,err)&
       result(o)
    DOUBLE_T, intent(in):: a, b, err!数値a,bと許容誤差err
    o = abs(a - b) < err
  end function nearlyEquals

!浮動小数点型複素数から、誤差を含めた比較を用いて整数型複素数を返す
  type(intcomp) function complexguess(z) result(iz)
#define RE real(z)
#define IM aimag(z)
#define GMAX 100
    DOUBLE_T, parameter:: err = 1e-8!許容誤差
    COMPLEX_T, intent(in):: z!浮動小数点型複素数
    integer:: inum, iden, gcd_
    DOUBLE_T:: dabs, darg
    dabs = abs(z)!複素数の実数ノルム
    darg = atan2(IM, RE)!複素数の実数位相(弧度法)
!整数ノルムの分子と分母をループさせることで実数のパターンをつくり、この実数ノルムと一致する整数ノルムを探す。分解能は1/GMAX
    ABSLOOP_NUM: do inum = -GMAX, GMAX!整数ノルムの分子のループ
       ABSLOOP_DEN: do iden = 1, GMAX!整数ノルムの分母のループ
          if(nearlyEquals(dabs, dble(inum)/iden, err)) then!type==1(|z|=num/den)に等しい場合
             gcd_ = getgcd(abs(inum), iden)
             iz%abs_num = inum/gcd_
             iz%abs_den = iden/gcd_
             iz%abs_type = 1
             go to 10
          else if (nearlyEquals(dabs, sqrt(inum+0d0)/iden, err)) then!type==2(|z|=sqrt(num)/den)に等しい場合
             gcd_ = getgcd(abs(inum), iden)
             iz%abs_num = inum/gcd_
             iz%abs_den = iden/gcd_
             iz%abs_type = 2
             go to 10
          end if
       end do ABSLOOP_DEN
    end do ABSLOOP_NUM
    iz%abs_type = 0
    return
!整数位相の分子と分母をループさせることで実数のパターンをつくり、実数位数と一致する整数位数を探す。分解能は1/GMAX
10  ARGLOOP_NUM: do inum = -GMAX, GMAX!整数位相の分子のループ
       ARGLOOP_DEN: do iden = 1, GMAX!整数位相の分母のループ
          if(nearlyEquals(darg, dble(inum)/iden*pi, err)) then!t = num/den*piのとき
             gcd_ = getgcd(abs(inum), iden)
             iz%arg_num = inum/gcd_
             iz%arg_den = iden/gcd_
             iz%arg_type = 1
             go to 20
          end if
       end do ARGLOOP_DEN
    end do ARGLOOP_NUM
    iz%arg_type = 0
    return
20  if(iz%arg_num < 0) then!num<0のとき、type==2(t=-num/den)としてnumの符号を反転する。
       iz%arg_num = -iz%arg_num
       iz%arg_type = 2
    end if
  contains
!2つの整数a,bの最大公約数を求める。
    integer function getgcd(a, b) result(c)
      integer, intent(in):: a, b
      integer:: n
      c = 1
      do n = 1, min(a,b)
         if(mod(a,n) == 0 .and. mod(b,n) == 0) c = n
      end do
    end function getgcd
#undef RE
#undef IM
  end function complexguess
end module m_commons
