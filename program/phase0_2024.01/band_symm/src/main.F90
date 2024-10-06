!TSPACEを用いたPHASEバンド簡約プログラムband_symm
!Copyright (C) 2014-2015 Ryosuke Tomita, Hiroyuki Oshima and Mineo Saito
!TSPACE: Copyright (C) 1995 Akira Yanagi (http://aquarius.mp.es.osaka-u.ac.jp/~tspace/)
!PHASE:  Copyright (C) 2002-2014 NIMS (https://azuma.nims.go.jp/software/phase)

#define READSTOP(func, fileno) \
if(.not.func(fileno)) stop 'error while '//"func";\
rewind(fileno)

program main
  use m_commons, only: stdout, openInputFromArgOrDefaultPath,openOutput
  use m_tspace_defines, only: character_table
  use m_tspace_commons, only: set_lattice_parameters, set_kpoint_irep, get_kgroup, set_pg_name,set_class,print_kpoint_and_kgroup
  use m_container_commons, only: readline_until_tag, k_state, c_state
  use m_container_lattice, only: trp2b, get_cstate, read_lattice_system, read_lattice_parameter, read_lattice_vector, read_kpoint_data
  use m_container_psicoef, only: ikpoint, nspin, readmagstate_psicoef, readstate_psicoef, ndim_spinor
  use m_classification, only: calcAll, printAll, deleteAll
  use m_compatibility, only: add_kpoint, print_compat_table, isDuplicateKpoint
  use m_character_table, only: print_character_table_for_TeX, get_character_table, print_character_table
  implicit none

!入出力用装置番号
  integer, parameter:: fileno = 1000
!出力用装置番号(reduce.data,character.tex)
  integer, parameter:: freduce(2) = (/1001,1011/), ftex = 1002
!k点の総数とスピンインデックス
  integer:: nikps, is, aa
!kpoint.dataから読み込んだk点のリストikps(1:4,1:nikps)
  integer, allocatable:: ikps(:,:)
!空間群の情報
  type(c_state), target:: cst
!k点及びk群の状態
  type(k_state), target:: kst
!指標表
  type(character_table), target:: cht
!**m.saito**
  WRITE(6,*) '** band_symm 202112'
  WRITE(6,*) '** MODIFICATION: main.F90, reduce_band.F90'
  WRITE(6,*) '   degeneracy.F90,tspace_commons.F90 band_symm.pl'

!インプットファイル1(kpoint.data)を開く
  if(.not.openInputFromArgOrDefaultPath(fileno, './kpoint.data'))&
       stop 'error while opening kpoint.data'
!kpoint.dataの読み込み、nikpsとikpsを決める
  call read_kpoint_data(fileno, nikps, ikps)
  close(fileno)

!インプットファイル2(band_sym_input.data)を開く
  if(.not.openInputFromArgOrDefaultPath(fileno, './band_sym_input.data'))&
       stop 'error while opening band_sym_input.data'

!格子情報、パラメータ、セルベクトルの読み込み
  READSTOP(read_lattice_system, fileno)
  READSTOP(read_lattice_parameter, fileno)
  READSTOP(read_lattice_vector, fileno)
  cst = get_cstate()
  call set_lattice_parameters(cst%latparam)

!空間群出力
  write(stdout, '(a)') '#START OF SPACE GROUP'
  call tspgds !TSPACE関数
  write(stdout, '(a)') '#END OF SPACE GROUP'
  write(stdout, '(a)') '#START OF K-POINTS'

!band_sym_input.dataの波動関数の情報が始まる行まで読み飛ばし
  if(readline_until_tag(fileno, '##PSIINPSTART') < 1 .or. &
       .not.readmagstate_psicoef(fileno))&
       stop 'read while reading PSIINPSTART'

!磁性を含めた計算であるかどうか
  kst%kpoint%isdouble = nspin > 1
!基本逆格子をブラベ逆格子に変換する行列(3x3)
  kst%trp2b = trp2b

!磁性を含めた計算のとき、reduce_up.dataとreduce_down.dataを出力用に開く。非磁性のとき、reduce.dataのみ出力用に開く。
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
 if(ndim_spinor == 2) then
     if(.not.openOutput(freduce(1), './reduce_soc.data'))&
          stop 'error while opening reduce.data'
  else if(nspin == 2) then
!--------------------------------------------------------------------------------
     if(.not.openOutput(freduce(1), './reduce_up.data'))&
          stop 'error while opening reduce.data'
     if(.not.openOutput(freduce(2), './reduce_down.data'))&
          stop 'error while opening reduce.data'
  else
     if(.not.openOutput(freduce(1), './reduce.data'))&
          stop 'error while opening reduce.data'
  end if

!指標表の出力ファイルを開く
  if(.not.openOutput(ftex, './character.tex'))&
       stop 'error while opening character.tex'

!TeXのプリアンブル文を記述
  write(ftex,*) '\documentclass{article}'
  write(ftex,*) '\begin{document}'
  write(ftex,*) '\section{Character table}'
!band_sym_input.dataの各k点の状態が記述されている行まで読み飛ばし
  do while(readstate_psicoef(fileno))
!k点および点群の情報を変数に格納
     if(ikpoint > nikps) cycle
     kst%kpoint%index = ikpoint
     kst%kpoint%comp = nint(matmul(dble(ikps(1:3, ikpoint)), kst%trp2b))
     kst%kpoint%common_denom = ikps(4,ikpoint)
!k点の情報からk群と既約表現を取得する
     call set_kpoint_irep(kst%kpoint,nspin)
     kst%kgroup = get_kgroup()!this function set kgroup order
     call set_pg_name(kst%kgroup)
     call set_class(kst%kpoint, kst%kgroup)!this function set class number of class
     !call tspgcr
!k点と群の名前の出力
     call print_kpoint_and_kgroup(stdout, kst%kpoint, kst%kgroup)

!指標表の取得
     cht = get_character_table(kst%kpoint, kst%kgroup)!this function set number of representation
     !call print_character_table(cht, kst%kgroup, stdout)
!k-pathの中で名前が重複しないk点のとき、適合関係を計算するk点のリストに登録し、指標表を出力する
     if(.not.isDuplicateKpoint(kst%kpoint)) then
        call add_kpoint(kst%kpoint, cht)
        call print_character_table_for_TeX(cht, kst%kgroup, kst%kpoint, ftex)
     end if

!波動関数の簡約計算を行い、スピンごとに出力する
     call calcAll(kst, cht, fileno)
     do is = 1, nspin/ndim_spinor
        call printAll(kst, is, freduce(is))
     end do
     call deleteAll
  end do
!k点の簡約計算の終わり
  write(stdout, '(a)') '#END OF K-POINTS'
  close(stdout)
  open(stdout, status='scratch')

!ファイルのクローズ及びデータ領域の開放
  if(allocated(ikps)) deallocate(ikps)
  close(fileno)
  do is = 1, nspin/ndim_spinor
     close(freduce(is))
  end do

!TeXファイルに適合関係を出力する
  write(ftex,*) '\clearpage'
  write(ftex,*) '\section{Compatibility table}'
  call print_compat_table(ftex)
  write(ftex,*) '\end{document}'
  close(ftex)
  close(stdout)
  stop
end program main
