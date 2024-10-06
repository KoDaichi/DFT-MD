!TSPACE関数のラッピング、以下は構造体や定数の定義部分
!Copyright (C) 2014-2015 Ryosuke Tomita 

module m_tspace_defines
#include "defines.F"
  implicit none
!結晶学で群論を扱う上での定数
  integer, parameter:: num_operator_max = 48, num_spacegroup_max = 230,&
!既約表現の数の最大値
       num_imr_max = 12, num_class_max = 12, &
!生成元の数の最大値
       num_generator_max = 3,&
!類の位数の最大値
       num_operator_per_class_max = 8,&
!表現行列の最大次数
       dimmax = 6

!TSPACE関数tskfbzの返り値、ISINTERNALであるとき、そのk点はBZ内部に、ISOUTERのとき外部に存在する
  integer, parameter:: ISINTERNAL = 1, ISOUTER = 0

!TSPACE内の対称操作の文字列表現
  character(len=4), parameter:: &
!反転操作
       element_name_inversion = 'I',&
!Oh群での対称操作
       element_name_cubic(num_operator_max/2) =&
       (/'E   ','C2X ','C2Y ','C2Z ',&
       'C31+','C32+','C33+','C34+',&
       'C31-','C32-','C33-','C34-',&
       'C2A ','C2B ','C2C ','C2D ','C2E ','C2F ',&
       'C4X+','C4Y+','C4Z+','C4X-','C4Y-','C4Z-'/),&
!D6h群での対称操作
       element_name_hexagonal(num_operator_max/4) =&
       (/'E   ','C6+ ','C3+ ','C2  ','C3- ','C6- ',&
       'C211','C221','C231','C212','C222','C232'/)

!格子系のパラメータ
  type lattice_parameter
!格子系[PFICRH]
     character(len=NameLen):: lattice_type
     integer::&
!生成元の数
          num_generator,&
!生成元の回転操作のインデックス
          index_rotation(num_generator_max),&
!生成元の並進ベクトル
          translation(2, 3, num_generator_max)
!反転対称をもつかどうか
     logical:: isinversion
!格子パラメータ
     DOUBLE_T:: a, b, c, alpha, beta, gamma
  end type lattice_parameter

!TSPACE上で扱われるk点の情報
  type kvector
     integer::&
!k点のインデックス
          index,&
!k点の成分表示ijk(b表示)
          comp(3),&
!k点の共通分母N, kx = i/N*b1x + j/N*b2x + k/N*b3x
          common_denom,&
!comm2bz: k点とBZとの関係、内部の点(ISINTERNAL)か外部(ISOUTER)か境界(>1)か
          comm2bz
!TSPACEで名付けられるk点特有の名前
     character(len=2):: name
!isdouble: 二重表現であれば真、でなくば偽
     logical:: isdouble
  end type kvector

!適合関係に関するTSPACEから得られる情報
  type comtable
     type(kvector):: kp(2)
!群(small)Gおよび部分群(large)Hの既約表現の数
     integer:: nimr(2),&
!Hの指標をGの指標で簡約したときの係数
          comp(num_imr_max, num_imr_max)
  end type comtable

!既約表現
  type irrep
!既約表現の名前(m_characte_tableで名付ける)
     character(len=NameLen):: name, name_mulliken
!次元
     integer:: dim
!既約表現の、対称操作における指標
     COMPLEX_T:: comp(num_class_max),&
!表現行列
          wd(dimmax, dimmax, num_operator_max)
  end type irrep

!TSPACEから得られる指標表の情報
  type character_table
!群の位数および既約表現の数
     integer:: nclass, numimr, numop, maxdim
!群に属する既約表現
     type(irrep):: imrs(num_imr_max)
     character(len=NameLen):: class_name(num_class_max)
  end type character_table

!対称操作
  type operat
!対称操作の名前(m_tspace_commons::get_element_nameで名付けられる)
     character(len=NameLen):: name
     DOUBLE_T::&
!回転行列（ブラベ格子表示）
          rot(nxyz_, nxyz_),&
!並進ベクトル（ブラベ格子表示）
          tr(nxyz_),&
          vaxis(nxyz_), vang
     COMPLEX_T:: spinrot(2,2)
     integer:: index_for_gamma, index_for_all
  end type operat
  
!類
  type class
!類の名前(m_tspace_commons::get_class_nameで名付けられる)
     character(len=NameLen):: name

     integer::&
!類の位数
          order,&
!類に属する対称操作のインデックス
          ind(num_operator_per_class_max)
  end type class

!TSPACEから得られる対称操作の情報
!
  type group
     integer::&
!order: 群の位数
          order,&
!nclass:類の数
          nclass,&
!群が反転操作を含むとき、その対称操作の番号を、含まないとき0以下の値をもつ
          inversion
!群に属する類
     type(class):: cls(num_class_max)
!群に属する対称操作
     type(operat):: op(num_operator_max)
     character(len=NameLen)::&
!点群の名前(Schönflies notation)
          schname,&
!点群の名前(Hermann–Mauguin notation)
          hmname
  end type group
end module m_tspace_defines
