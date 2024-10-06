!PHASEから入力された空間群の生成元および格子の情報を読み込み、保存するモジュール
!Copyright (C) 2014-2015 Ryosuke Tomita

module m_container_lattice
#include "defines.F"
  use m_commons, only: stdout, unit_matrix, pi, pi2, euler2vector, get_recip_vector
  use m_tspace_defines, only: kvector, group, lattice_parameter, num_generator_max
  use m_tspace_commons, only: get_group, get_kgroup, get_kgroup_from_gammagroup, set_kpoint_irep, set_pg_name, set_class
  use m_container_commons, only: c_state, k_state, linebuf, parser_tag, readline, readline_until_tag
  
  implicit none
!m_container_common::readline_until_tag用のタグ
  character(*), parameter::&
       kvector_tag = 'KPOINT', &
       lattice_system_tag = 'LATTICE_SYSTEM',&
       generators_tag = 'GENERATORS',&
       lattice_parameter_tag = 'BRAVAIS_LATTICE_PARAMETER',&
       lattice_vector_tag = 'PRIMITIVE_LATTICE_PARAMETER'
  
!k点を格納する変数
  type(kvector), private, target:: kpoint
!格子定数の情報を格納する
  type(lattice_parameter), private, target:: latparam
!k点群
  type(group), private, target:: kgroup
!格子変換行列、実格子または逆格子のブラベ格子及びプリミティブ格子
  DOUBLE_T, dimension(3,3), private, target:: alat_pri, alat_bra, blat_pri, blat_bra       
  DOUBLE_T, dimension(3,3), public, target:: trp2b = unit_matrix, &
       trb2p = unit_matrix, trp2b_a = unit_matrix, trb2p_a = unit_matrix
contains

!読み込んだ情報およびTSPACEの内部変数をもとにk点およびk点群の情報をモジュール内の変数に格納する
  subroutine set_kpoint_and_kgroup()
    logical, save:: isgammaset = .false.
    type(kvector), save:: k_gamma!GAMMA点
    type(group), save:: gr_gamma!GAMMA点の群
!GAMMA点の群を取得する
    if(.not.isgammaset) then
       k_gamma%comp = 0
       k_gamma%common_denom = 1
       call set_kpoint_irep(k_gamma)
       gr_gamma = get_kgroup()
       isgammaset = .true.
    end if
!k点の成分をブラベ逆格子のものに変換
    kpoint%comp = nint(matmul(kpoint%comp, trp2b))
    call set_kpoint_irep(kpoint)
!GAMMA点の群から必要な対称操作だけを取得して、k群を用意する
    kgroup = get_kgroup_from_gammagroup(gr_gamma)
    call set_pg_name(kgroup)
    call set_class(kpoint, kgroup)
  end subroutine set_kpoint_and_kgroup

!kpoint.dataの読み込み.k点の総数をnkpに、k点の成分をkps(1:4,1:nkp)に記憶する
!kps(1:4,i)はi番目のk点の成分分子x,y,zと共通分母Nを記憶する。
  subroutine read_kpoint_data(nfin, nkp, kps)
    integer, intent(in):: nfin
    integer, intent(out):: nkp
    integer, allocatable, intent(out):: kps(:,:)
    integer:: ikp
    if(readline(nfin) < 1) then
       nkp = 0
       return
    endif
    read(linebuf, *) nkp
    allocate(kps(4,nkp))
    ikp = 0
    do while(readline(nfin) > 0 .and. ikp < nkp)
       read(linebuf, *) kps(1:4, ikp+1)
       ikp = ikp + 1
    enddo
  end subroutine read_kpoint_data

!k点の読み込み
  logical function read_kpoint(nfin, reset, notag)
    integer, intent(in):: nfin
    logical, intent(in), optional:: reset, notag
    integer:: ik
    integer, save:: count = 1
    read_kpoint = .false.
    
    if(present(reset) .and. reset) count = 1
!読み込み開始処理、m_container_psicoef::read_psicoefと同様
!notagタグがないときに、KVECTORタグを読み込むまでファイル内容を読み飛ばす
    if(.not.(present(notag) .and. notag) .and. &
         readline_until_tag(nfin, kvector_tag) < 1) return
!読み込み終了処理、m_container_psicoef::read_psicoefと同様
    do ik = 1, count
       if(READERR(nfin, parser_tag)) then
          count = 1
          return
       end if
       read(linebuf, *) kpoint%index, kpoint%comp, kpoint%common_denom, kpoint%isdouble
       kpoint%name = ""
    enddo
    count = count + 1
    read_kpoint = .true.
  end function read_kpoint
  
!結晶系および生成元、反転対称性の読み込み
  logical function read_lattice_system(nfin)
    integer, intent(in):: nfin
    integer:: igen

    read_lattice_system = .false.
!読み込み開始処理、m_container_psicoef::read_psicoefと同様
    if(readline_until_tag(nfin, lattice_system_tag) < 1) return
    if(READERR(nfin, parser_tag)) return
    read(linebuf, *) latparam%lattice_type, latparam%isinversion

    if(readline_until_tag(nfin, generators_tag) < 1) return
!生成元の読み込み、最大3つ
    do igen = 1, num_generator_max
       if(READERR(nfin, parser_tag)) exit
       read(linebuf, *) latparam%index_rotation(igen), latparam%translation(:, :, igen)
       latparam%num_generator = igen
    enddo
    read_lattice_system = .true.
  end function read_lattice_system

!格子パラメータの読み込み
  logical function read_lattice_parameter(nfin)
    integer, intent(in):: nfin

    read_lattice_parameter = .false.
!読み込み開始処理、m_container_psicoef::read_psicoefと同様
    if(readline_until_tag(nfin, lattice_parameter_tag) < 1) return

!読み込み終了処理、m_container_psicoef::read_psicoefと同様
    if(READERR(nfin, parser_tag)) return
    read(linebuf, *) latparam%a, latparam%b, latparam%c
    if(READERR(nfin, parser_tag)) return
    read(linebuf, *) latparam%alpha, latparam%beta, latparam%gamma
    read_lattice_parameter = .true.
  end function read_lattice_parameter

!プリミティブ格子の読み込み(データの順番は1カラムから逆格子、4カラムから実格子)
  logical function read_lattice_primitive(nfin)
#define READLATTICE(num)\
    read(linebuf, *) blat_pri(num, :), alat_pri(num, :)

    integer, intent(in):: nfin
    read_lattice_primitive = .false.
    if(readline_until_tag(nfin, lattice_vector_tag) < 1)&
         return
    if(READERR(nfin, parser_tag)) return
    READLATTICE(1)!a_1の読み込み
    if(READERR(nfin, parser_tag)) return
    READLATTICE(2)!a_2の読み込み
    if(READERR(nfin, parser_tag)) return
    READLATTICE(3)!a_3の読み込み
    read_lattice_primitive = .true.
#undef READLATTICE
  end function read_lattice_primitive

!TSPACEの内部変数からブラベ格子ベクトルを設定する
  subroutine set_lattice_braveis
    DOUBLE_T, save:: t0 = pi/180d0
    alat_bra = euler2vector(latparam%alpha * t0,&
         latparam%beta * t0,&
         latparam%gamma * t0&
         )
    alat_bra(:, 1) = alat_bra(:, 1) * latparam%a
    alat_bra(:, 2) = alat_bra(:, 2) * latparam%b
    alat_bra(:, 3) = alat_bra(:, 3) * latparam%c
    blat_bra = get_recip_vector(alat_bra)
  end subroutine set_lattice_braveis

!読み込んだプリミティブ格子及びブラベ格子から、ブラベ格子表示のベクトルとプリミティブ格子表示のベクトルを相互変換する格子変換行列(trp2b, trb2p)を用意する
  subroutine set_trp2b
    trb2p = matmul(blat_bra, transpose(alat_pri)) /pi2! blat_bra * blat_pri^-1
    trp2b = matmul(blat_pri, transpose(alat_bra)) /pi2! blat_pri * blat_bra^-1
    trb2p_a = transpose(trp2b)
    trp2b_a = transpose(trb2p)
  end subroutine set_trp2b
 
!プリミティブ格子の読み込み、ブラベ格子および格子変換行列の設定
  logical function read_lattice_vector(nfin)
    integer, intent(in):: nfin
    read_lattice_vector = read_lattice_primitive(nfin)
    if(.not.read_lattice_vector) return
    call set_lattice_braveis
    call set_trp2b
  end function read_lattice_vector

  subroutine trans_transvector
    integer:: i, j
    DOUBLE_T:: v(3)
#ifdef DEBUG
    do i = 1, latparam%num_generator
       write(stdout, *) (latparam%translation(:, j, i), j = 1,3)
    end do
#endif
    do i = 1, latparam%num_generator
       v(:) = dble(latparam%translation(1,:,i))
       latparam%translation(1,:,i) = nint(matmul(v, trp2b_a))
       do j = 1, 3
          latparam%translation(1,j,i) = mymod(latparam%translation(1,j,i),&
               latparam%translation(2,j,i))
       end do
    end do
#ifdef DEBUG
    do i = 1, latparam%num_generator
       write(stdout, *) (latparam%translation(:, j, i), j = 1,3)
    end do
#endif
  contains
    integer function mymod(a,b)
      integer, intent(in):: a,b
      mymod = mod(abs(a), b)
      if(a < 0) mymod = -mymod
    end function mymod
  end subroutine trans_transvector

!得られたk点の情報をコピーし、それを返す
  type(k_state) function get_kstate() result(ks)
    ks%kpoint = kpoint
    ks%kgroup = kgroup
    ks%trp2b = trp2b
    ks%trb2p = trb2p
  end function get_kstate

!得られた結晶全体の情報をコピーし、それを返す
  type(c_state) function get_cstate() result(cs)
!    call trans_transvector
    cs%latparam = latparam
    cs%crystal_group = get_group()
    cs%trp2b = trp2b_a
    cs%trb2p = trb2p_a
  end function get_cstate
end module m_container_lattice
