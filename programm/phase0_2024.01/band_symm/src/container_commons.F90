!ファイルのリーディングおよび計算で使用する変数の格納に関するモジュールのひとつ、ファイルリーディングに関する共通の変数および関数を定義する
!Copyright (C) 2014-2015 Ryosuke Tomita

#define READSELECT(nfin) \
  select case(readline(nfin));\
  case(-1); return;\
  case(0); exit;\
  end select

  module m_container_commons
    use m_tspace_defines, only: kvector, group, lattice_parameter
#include "defines.F"
    implicit none
!ファイルの読み込みバッファ
    character(len = 255), public:: linebuf
!ファイルのデータ区切り(parser)およびコメント文(comment)の始まりを表す文字列
    character(*), parameter::&
         parser_tag = '#', &
         comment_tag = '!'
!一つの縮退点に入るバンドの総数
!SAITO
!    integer, public, parameter:: max_num_deg_band = 3
    integer, public, parameter:: max_num_deg_band = 4
    type energy_tag
       integer:: ie
!エネルギー準位の、基準からの数値(UP,DOWN)
       DOUBLE_T:: hartree
    end type energy_tag
    
    type degeneracy
!一つの縮退点における準位の数およびそのインデクス
       integer:: ne, ie(max_num_deg_band)
    end type degeneracy
    
    type gvector_t
       integer:: num
       integer, pointer:: comp(:, :)
    end type gvector_t
    
    type psicoef_t
       integer:: num_energy, nspin, num_deg_point(2)!, nspinor
       type(gvector_t):: gv
       type(energy_tag), pointer:: e_tag(:,:)
       type(degeneracy), pointer:: deg(:,:)
       COMPLEX_T, pointer:: comp(:,:,:)
    end type psicoef_t
    
    type k_state
       type(kvector):: kpoint
       type(group):: kgroup
       DOUBLE_T:: trp2b(3,3), trb2p(3,3)
    end type k_state
    
    type c_state
       type(lattice_parameter):: latparam
       type(group):: crystal_group
       DOUBLE_T:: trp2b(3,3), trb2p(3,3)
    end type c_state
    
    type transd_vector_index
       integer:: numop, numg
       integer, pointer:: gv(:,:), kv(:)
    end type transd_vector_index
    
  contains
!ファイルnfinから意味のある行(コメントタグから始まらない、かつ空でない)を一行読み込み、バッファに格納する関数、リーディングエラーで-1,EOFで0、正常終了で1を返す
    integer function readline(nfin)
      integer, intent(in):: nfin
      do while(.true.)
         read(nfin, fmt = '(a)', err = 101, end = 100) linebuf
         linebuf = trim(linebuf)
         if(len(linebuf) > 0 .and. linebuf(1:1) /= comment_tag) exit
      enddo
      readline = 1; return
101   readline = -1; return
100   readline = 0; return
    end function readline

!指定されたtagがある行を読み込むまでファイルnfinの行を読み飛ばす。返り値の意味はreadlineのものと同じ。正常終了の時、バッファにはtagがある行が格納される
  integer function readline_until_tag(nfin, tag) result(iout)
    integer, intent(in):: nfin
    character(*), intent(in):: tag

    iout = 1
    do while(.true.)
       select case(readline(nfin))
       case(-1)
          iout = -1; return
       case(0)
          iout = 0; return
       end select
       if(index(linebuf, tag) > 0) exit
    end do
    
  end function readline_until_tag
end module m_container_commons
