#!/usr/bin/perl
#reduce.dataから既約表現付きのバンド図を生成するband_symm.pl
#Powerer by Ryosuke Tomita at 2015/02/11

use strict;
use warnings;
use Getopt::Long;

#グローバル変数
#環境に応じて固定する変数
my $HAR2EV = 27.2113834;#ハートリー単位とeVの変換、1hartree = 27.2112eV
my $GNUPATH = "gnuplot";#gnuplotのパス
my %TSPACE_KNAME_TABLE =
	("GM" => "{/Symbol G}",
	 "DT" => "{/Symbol D}",
	 "SM" => "{/Symbol S}",
	 "LD" => "{/Symbol L}");#TSPACEで出力されるk点の名前(key)とepsで用いるギリシャ文字(value)との対応

#レイアウト用設定(一部オプションから指定できる)
my @IMRFONT = ("Helvetica,12", "Helvetica,12");#既約表現ラベルのデフォルトフォント、アップスピン(or 非磁性軌道)とダウンスピン軌道の2成分
my $TICSFONT = "Helvetica,12";#目盛ラベルのデフォルトフォント
my @IMRCOLOR = ('black', 'gray');#既約表現の文字色、アップスピン(or 非磁性軌道)とダウンスピン軌道の2成分
my @LINECOLOR = ('black', 'gray');#バンドの線の色、アップスピン(or 非磁性軌道)とダウンスピン軌道の2成分
my $YLABEL = 'Energy (eV)';#y軸のラベル
my $PLOTTYPE = "lines";# バンドをプロットするときの、gnuplotのオプション{solid_circles|lines}
my %MARGINSET = ('t' => 3, 'b' => 3);#マージンを変更する箇所とその値。topとbottomだけ変更し、leftとrightはgnuplotのデフォルトのマージンを使用する
my $ENERGY_IMR_DUPLICATE_THRESHOLD = 0.0001;#偶然縮退すると判断する二つのバンドにおけるエネルギー差の閾値(eV)
my @nonconventional= (); # NONCONVENTIONAL IMR =1


#隠しオプション
my $IS_ENUM_IMR = 1; #既約表現が重なるとき、横に並べて配置するかどうか
my $OFFSET = 1; #バンド図両端の既約表現は縦軸と重なり見難い。その既約表現を少し中央寄りにどれだけシフトするか。デフォルトは1k点分





#各種デフォルト設定(一部オプションから変更できる)
my $IS_DROW_IMR = 1;#既約表現のラベルを記述するかどうか(1:True,0:False)
my $IS_DROW_FERMILINE = 0;#フェルミ準位(y=0)に破線を引くかどうか(1:True,0:False)
my $IMR_TYPE = 1;#既約表現の種類(0:TSPACE表現,1:Mulliken表記)
my $SW_BAND = 1;#バンド交差機能を実行するかどうか(1:True,0:False)
my $IS_UNLINK_SRC = 1;#FILE_BANDとFILE_GPIを実行後に削除するかどうか(1:True,0:False)
my $K_GROUP_TYPE = 1;#k点の群の名前の種類(1のときシェーンフリース、2のとき国際記号を表示)
my $IS_KNAME_EDGE_DOWN = 1;#k-pathの節と節間のk点の名前を字下げするかどうか(1:True,0:False)
my $EINC_DEVIDE_NUM = 10;#EINCを(EMAX-EMIN)/Nで表現するときの整数N、つまりy軸の大目盛りの数
my $ENERGY_CROSS_PRECISION = 0.01;#接近していると判断する2バンドのエネルギー差の閾値(eV)2バンドのエネルギー差がこの値未満になるとき、この2バンドが交差しているとみなす。
my $NUMIMR = 1.0;#k-pathの節と節の間にいくつ規約表現をつけるか。 
#出力ファイルパス
my $OUTFILE = "./band_reduce.eps";#出力するバンド図のデフォルトパス
my $FILE_BAND = "./$$.dat";#固有値を出力するファイルのパス、デフォルトで実行後に削除される
my $FILE_GPI = "./$$.gpi";#gnuplot用のスクリプトを出力するパス、デフォルトで実行後に削除される

#入力データ群
my @ENERGY = ();#固有値
my @bufE = ();  #ENERGYのバッファデータ
my @DIME = (); #固有値の縮退度
my ($NKVEC, $NBAND, $NSPIN, $EFERMI) = (0,0,1,0.0);#k点の数、バンドの数、スピンの数、フェルミ準位(eV単位)
my @KVEC = ();#各k点の座標(整数座標、分子x,y,zが要素0,1,2に、共通分母Nが要素3に保存される。)

my @IMRNAME = ();#既約表現の名前
my @IMRINDEX = ();#既約表現の番号
my @K_NAME = ();#k点の名前
my @K_GROUP_NAME = ();#k点の群の名前
my @MAXDIM = ();#各k点のバンドにおける縮退数の最大値

my @DK = ();#k-pathの重み付き道のり
my @LK = ();#既約表現を表示しないk点の場合は0、表示するk点でかつk-pathの節である場合は1、節でないときは2
my @NUM_BAND_K = ();#各k点で、$EMINと$EMAXで区切られたエネルギー範囲に表示するバンドの数
my ($EMIN,$EMAX,$EINC) = (0,0,0);#バンドを表示するエネルギー値の最小最大値と目盛り(eV)
my @IMRFONTVALUE = (12,12);
my %OPTS = ();#オプション

print "** band_symm.pl 2022 (Ver.3) ** \n";
#オプション解析
GetOptions(\%OPTS, 
	'erange=s', #表示するエネルギーレンジ
	'einc=s', #エネルギー目盛り
	'with_fermi=s',#フェルミ準位の指定
	'linecolor=s',#バンドのカラー表示
	'imrcolor=s',#既約表現の文字の色
	'imrfont=s',#既約表現ラベルのフォント
	'numimr=s',#節-節間に表示する規約表現の数(デフォルト1)
	'ticsfont=s',#軸ラベルのフォント
	'imrtype=s',#既約表現の記号の種類
	'nonecross',#バンド交差機能の無効化
	'kgrouptype=s',#k点の群の名前の種類
	'offset=s',	   
	'h',#ヘルプ
);
if($OPTS{'h'} || @ARGV<1) {#ヘルプオプションを設定したか、reduce.dataが指定されていないときヘルプを表示して終了
    &help;
    exit;
}


if($OPTS{'imrtype'} && ($OPTS{'imrtype'} =~ /^[nN]/)){
	$IMR_TYPE = 0;
}
if($OPTS{'kgrouptype'} && ($OPTS{'kgrouptype'} =~ /^[hH]/)){
	$K_GROUP_TYPE = 2;
}

$SW_BAND = !$OPTS{'nonecross'};

if($OPTS{"imrfont"}){
	my @tics = split(/,/, $OPTS{'imrfont'});
	for(my $i = 0; $i < @tics; $i+=2){
		$IMRFONT[$i/2] = join(',',$tics[$i],$tics[$i+1]);
		$IMRFONTVALUE[$i/2] = $tics[$i+1];
	}
}
if($OPTS{"offset"}){
	$OFFSET = $OPTS{'offset'};
}
if($OPTS{"numimr"}){
	$NUMIMR = $OPTS{'numimr'};
}
if($OPTS{"ticsfont"}){
	$TICSFONT = $OPTS{'ticsfont'};
}
if($OPTS{'erange'}){
	($EMIN,$EMAX) = split(/\s*,\s*/, $OPTS{'erange'});
}
if($OPTS{'einc'}){
	$EINC = $OPTS{'einc'};
}
if($OPTS{'with_fermi'}){
	$EFERMI = &parseFermi($OPTS{'with_fermi'});
	$IS_DROW_FERMILINE = 1;
}
if($OPTS{'imrcolor'}){
	my @color = split(/,/,$OPTS{'imrcolor'});
	for(my $i = 0; $i < @color; $i++){
		$IMRCOLOR[$i] = $color[$i];
	}
}
if($OPTS{'linecolor'}){
	my @color = split(/,/,$OPTS{'linecolor'});
	for(my $i = 0; $i < @color; $i++){
		$LINECOLOR[$i] = $color[$i];
	}
}


#reduce.dataの読み取り

my @file_reduce = @ARGV;#$file_reduce[0]にアップスピン、$file_reduce[1]にダウンスピンのreduce.dataを指定
$NSPIN = @file_reduce;
for(my $i = 0; $i < $NSPIN; $i++){
	$IMRNAME[$i] = [];
	$IMRINDEX[$i] = [];
	open(my $IN, $file_reduce[$i]) or die $file_reduce[$i] ." is not found or busy\n";#ファイルが開くことができない場合は終了
	&readReduceData($IN, $i);#それぞれのreduce.dataを読み取り、各状態の固有値と既約表現を取得する
	close($IN);
}
$NKVEC = @KVEC;
$NBAND = @{$ENERGY[0][0]};
&makeDKandLK;#@DKと@LKを指定する
&band2ev;#固有値の単位と基準を変換
if($SW_BAND){
	&BandCrossSet;#バンドの交差を表現する
}

#バンド図作成スクリプトおよびデータの作成
my $OUT;
open($OUT, ">", $FILE_BAND) or die "$FILE_BAND is busy";
&getBandValues($OUT);#固有値データを作成し、ファイルに出力
close($OUT);

open($OUT, ">", $FILE_GPI) or die "$FILE_GPI is busy";
&getGnuScript($OUT);#Gnuplotスクリプトを作成し、ファイルに出力
close($OUT);

print `$GNUPATH $FILE_GPI 2>&1`;#Gnuplotの実行、コンソールからの出力をパイプする
if($IS_UNLINK_SRC){
	unlink $FILE_GPI, $FILE_BAND;
}else{
	print "SRCFILES:%s,%s\n",$FILE_GPI, $FILE_BAND;
}





######################################################################
#固有値をフェルミ準位基準にし、単位をeVに変換
######################################################################
sub band2ev{
	my $_min = 1000;
	my $_max = -1000;
	for(my $is = 0; $is < $NSPIN; $is++){
		for(my $ik = 0; $ik < $NKVEC; $ik++){
			for(my $ib = 0; $ib < $NBAND; $ib++){
				$ENERGY[$is][$ik][$ib] *= $HAR2EV;
				$ENERGY[$is][$ik][$ib] -= $EFERMI;
				$bufE[$is][$ik][$ib] = $ENERGY[$is][$ik][$ib];
				if($ENERGY[$is][$ik][$ib] < $_min){$_min = $ENERGY[$is][$ik][$ib];}
				if($ENERGY[$is][$ik][$ib] > $_max){$_max = $ENERGY[$is][$ik][$ib];}
			}
		}
	}
	unless($OPTS{'erange'}){#エネルギー表示オプションがないときはEMIN,EMAXをデータから判断する
		$EMIN = $_min-1; $EMAX = $_max+1;
	}
	unless($OPTS{'einc'}){#エネルギー目盛りオプションがないとき、EMIN,EMAXから決定する
		$EINC = int(($EMAX - $EMIN)/$EINC_DEVIDE_NUM);#y軸をEINC_DEVIDE_NUM分割する
		if($EINC == 0){
			$EINC = 0.1;
		}
	}
}





######################################################################
#バンドプロットのための固有値を出力する
######################################################################
sub getBandValues{
	my $OUT = shift;
	for(my $ik = 0; $ik < $NKVEC; $ik++){
		print $OUT $DK[$ik], " ";#DKを出力
		for(my $is = 0; $is < $NSPIN; $is++){
			print $OUT join(" ",@{$ENERGY[$is][$ik]})," ";#バンド固有値を出力
		}
		print $OUT "\n";
	}
}





######################################################################
#バンドプロットのためのgnuplotスクリプトを出力
######################################################################
sub getGnuScript{
	my $OUT = shift;
#プロット用の変数の定義
	print $OUT "TICSFONT = \"${TICSFONT}\"\n";
	for(my $i = 0; $i < $NSPIN; $i++){#スピン軌道ごとにフォントと色を分ける
		print $OUT "IMRFONT${i} = \"$IMRFONT[$i]\"\n";
		print $OUT "IMRCOLOR${i} = \"$IMRCOLOR[$i]\"\n";
		print $OUT "LINECOLOR${i} = \"$LINECOLOR[$i]\"\n";
	}

#軸ラベルおよびプロット範囲の設定
	print $OUT "set ylabel \"${YLABEL}\"\n";
	print $OUT "set yrange [$EMIN:$EMAX]\n";
	print $OUT "set xtics scale 0\n";
	print $OUT "set xtics font TICSFONT\n";
	print $OUT "set x2tics scale 0\n";
	print $OUT "set x2tics font TICSFONT\n";
	print $OUT "set ytics $EINC\n";
	foreach my $type (keys %MARGINSET){
	    print $OUT "set ${type}margin $MARGINSET{$type}\n";
	}
	
	my @tics = ([],[]);
	my $iik = 0;
	for(my $ik = 0; $ik < $NKVEC; $ik++){
		next unless($LK[0][$ik][1] > 0);
		push(@{$tics[0]}, sprintf("\"%s\" %lf",
			(($IS_KNAME_EDGE_DOWN && ($LK[0][$ik][1] != 1)) ? "\\n" : "").$K_NAME[$ik],  $DK[$ik]));
		push(@{$tics[1]}, sprintf("\"%s\" %lf",$K_GROUP_NAME[$ik],  $DK[$ik]));
	}
	print $OUT "set xtics (",join(",",@{$tics[0]}),")\n";
	print $OUT "set x2tics (",join(",",@{$tics[1]}),")\n";
	for(my $ik = 1; $ik < ($NKVEC-1); $ik++){
		next unless($LK[0][$ik][0] == 1);#kがk-pathの節でないときスキップ
		print $OUT "set arrow from ($DK[$ik]),($EMIN) to ($DK[$ik]),($EMAX) nohead lt 0\n";#k-pathの節で、バンド図に縦に線を引く
	}
	if($IS_DROW_FERMILINE){
		print $OUT "set arrow from 0.0,0.0 to 1.0,0.0 nohead lt 0\n";#フェルミ準位で破線を引く
	}

#既約表現ラベルの設定

	my $imrformat = "set label \"%s\" at first (%lf),(%lf) center font IMRFONT%d textcolor rgb IMRCOLOR%d\n";#既約表現ラベル宣言のフォーマット
	for(my $is = 0; $is < $NSPIN; $is++){
		my $imr = $IMRNAME[$is];
		next unless($IS_DROW_IMR);
		for(my $ik = 0; $ik < $NKVEC; $ik++){
			my $offfset =0;
			if( $ik ==0 ){
				$offfset = ($DK[1]-$DK[0])*$OFFSET ;
			}
			if($ik ==$NKVEC-1){
				$offfset = ($DK[$NKVEC-2]-$DK[$NKVEC-1])*$OFFSET;
				}
			next unless($LK[$is][$ik][0] > 0);#LKに登録されたikにだけラベルをつける。

			my @imrnew = ();#偶然縮退を考慮した新しい既約表現ラベル
			my @Isimr = ();
#2021/12 modify SAito
			my @posi_energy=(); #エネルギー位置を変化させる
#           my @count_degeneracy=();
            my $ENERGY_DIFFERENCE =0.0022*$IMRFONTVALUE[$is]*($EMAX-$EMIN);
#            
			for(my $ib = 0;$ib<$NBAND;$ib++){
			 $Isimr[$ib] = 0;
			}
#			print "$ENERGY_DIFFERENCE, $EMAX,$EMIN \n";
			for(my $ib = 0; $ib < $NBAND-1; $ib++){#偶然縮退を考慮して新しく既約表現のラベルを組みなおすためのバンドのループ
				#バンドibの固有値、ib+1の固有値
				my @en = ($bufE[$is][$ik][$ib], $bufE[$is][$ik][$ib+1]);
				#バンドibは表示可能かどうか、バンドib+1は表示可能かどうか
				my @esv = ($$imr[$ik][$ib] && ($en[0] <= $EMAX) && ($en[0] >= $EMIN),
						   $$imr[$ik][$ib+1] && ($en[1] <= $EMAX) && ($en[1] >= $EMIN));
				if(!$Isimr[$ib]){
					$imrnew[$ib]=$$imr[$ik][$ib];
				}
				$Isimr[$ib] ++;
				if($esv[0]){
					if($IS_ENUM_IMR){
					if($esv[1]){#ibが表示可能かつib+1も表示可能
				    	if(abs($bufE[$is][$ik][$ib] - $bufE[$is][$ik][$ib+1] ) < $ENERGY_IMR_DUPLICATE_THRESHOLD){
#						if(abs($bufE[$is][$ik][$ib] - $bufE[$is][$ik][$ib+1] ) <$ENERGY_IMR_DUPLICATE_THRESHOLD){
							#既約表現が接近するとき、二つの既約表現のx座標をずらす機能
							
							if($$imr[$ik][$ib] eq $$imr[$ik][$ib+1]){
#								if(abs($en[0] - $en[1]) <$ENERGY_IMR_DUPLICATE_THRESHOLD){
									$imrnew[$ib+1]=$imrnew[$ib];
									$imrnew[$ib]='';
#								}
							}else{
								if(!$Isimr[$ib+1]){
									$Isimr[$ib+1] = 1;
									$imrnew[$ib+1] =  $$imr[$ik][$ib+1];
									if($Isimr[$ib]==1){
										$imrnew[$ib] = sprintf("%s    ", $imrnew[$ib]);
										$imrnew[$ib+1] =sprintf("    %s", $imrnew[$ib+1]);
									}
									if($Isimr[$ib]==2){
										$imrnew[$ib+1] =sprintf("%s", $imrnew[$ib+1]);
										$Isimr[$ib+1]--;
									}
								}
							}
						}else{
							$imrnew[$ib+1] = $$imr[$ik][$ib+1];
						}
					}else{#ibが表示可能だがib+1は表示不可
						$imrnew[$ib+1] = '';
					}
					}
				}else{
					if($esv[1]){#ibは表示不可能だがib+1は表示可能
						$imrnew[$ib] = '';
						$imrnew[$ib+1] = $$imr[$ik][$ib];;
					}else{#ibとib+1どちらも表示不可能のとき
						$imrnew[$ib] = $imrnew[$ib] = '';
					}
				}
			}
#2021/12 modify Saito
#           3,4個のIMRが同じエネルギーを持つ場合
            for(my $ib = 0; $ib < $NBAND; $ib++){$posi_energy[$ib]=$bufE[$is][$ik][$ib];}
            for(my $ib = 2; $ib < $NBAND; $ib++)
            {  my $count_degeneracy=0;
               for (my $ibbb=0; $ibbb < $ib; $ibbb++)
               {  if($posi_energy[$ib] - $posi_energy[$ibbb] < $ENERGY_IMR_DUPLICATE_THRESHOLD and $imrnew[$ibbb] ne '')
                  {$count_degeneracy++;};
               }  
               if( $count_degeneracy > 1){$imrnew[$ib] =sprintf("              %s", $imrnew[$ib]);};
             }
#二つのレベルが縮退ではないが、接近しているときに、縦軸の値をずらして重ならないようにする         
            
            for(my $ib = 0; $ib < $NBAND-1; $ib++){	
			    $posi_energy[$ib+1]=  $bufE[$is][$ik][$ib+1];
			    if(abs($bufE[$is][$ik][$ib+1] - $bufE[$is][$ik][$ib]) < $ENERGY_IMR_DUPLICATE_THRESHOLD){  
                    $posi_energy[$ib+1]=$posi_energy[$ib]
                 ;}
                 else{if(($bufE[$is][$ik][$ib+1] - $posi_energy[$ib]) < $ENERGY_DIFFERENCE){
                     $posi_energy[$ib+1]=$posi_energy[$ib]+$ENERGY_DIFFERENCE
                 ;}
                ;}
            ;}       

			for(my $ib = 0; $ib < $NBAND; $ib++){$bufE[$is][$ik][$ib]=$posi_energy[$ib];
			                                     if($bufE[$is][$ik][$ib] > $EMAX-$ENERGY_DIFFERENCE/2){$imrnew[$ib]='';}
			                                     }
#END modify Saito
			
			for(my $ib = 0; $ib < $NBAND; $ib++){#組みなおした既約表現ラベルを出力するためのループ
				if($imrnew[$ib]){printf $OUT $imrformat, $imrnew[$ib], $DK[$ik]+$offfset, $bufE[$is][$ik][$ib], $is, $is;}
			}
#energy shift			
#		    	* appears above the baseline for NONCONVENTIONAL IMR
			my $ib=1;
			$imrnew[$ib]='*';
			$bufE[$is][$ik][$ib]=$EMIN + $ENERGY_DIFFERENCE;
			if($nonconventional[$ik] == 1){printf $OUT $imrformat, $imrnew[$ib], $DK[$ik]+$offfset, $bufE[$is][$ik][$ib], $is, $is;}
#		        * appears above the baseline for NONCONVENTIONAL IMR
		}
	}

#プロット部分の設定
	print $OUT "set term postscript eps enhanced color solid\n";
	print $OUT "set output \"$OUTFILE\"\n";
	print $OUT "set nokey\n";

	print $OUT "p ";
	for(my $is = 0; $is < $NSPIN; $is++){
		for(my $ib = 0; $ib < $NBAND; $ib++){
			printf $OUT "\"$FILE_BAND\" u 1:%d w $PLOTTYPE linewidth 1 lc rgb LINECOLOR%d,\\\n", $NBAND*$is+$ib+2, $is;
		}
	}
	
	print $OUT "NaN\n";
}





######################################################################
#reduce.dataの読み込み、データを$REDUCEと$ENERGY、$KVECに保存する
######################################################################
sub readReduceData{
	my ($IN,$is) = @_;
	my ($ib, $ik, $iimr) = (0,0,0);

	my $imr = $IMRNAME[$is];#既約表現の名前(二重配列。ik行ib列にik番目のk点とib番目のバンドにおける既約表現の名前が入る)
	my $imrnum = $IMRINDEX[$is];#既約表現の番号(二重配列。ik行ib列にik番目のk点とib番目のバンドにおける既約表現の番号が入る)バンド交差機能において、既約表現の比較に用いる
#簡約ファイル読み込み

	foreach my $line (<$IN>){
		$line =~ s/^\s*(.*)\s*\n?$/$1/;
#コメントまたは空白行のスキップ
		next if(!$line || ($line =~ /^\s*\#/));
#プロパティ行(!ではじまる)の読み込み
		if($line =~ /^!/){
			my @tag = ($line =~ /[^!\s]+/g);
			if(@tag < 4){#バンドインデックスの行
				$ib = $tag[0] -1;#インデックス
				my $de = $tag[1];#固有値(hatree)
				my $nd = $tag[2];#縮退度
				for(my $id = 0; $id < $nd; $id++){
					$ENERGY[$is][$ik][$ib+$id] = $de;
					$DIME[$is][$ik][$ib+$id] = $nd;
				}
				$iimr = 0;#既約表現のインデックスをリセット
			}else{#k点の成分、群の行
				$ik = $tag[0] -1;#インデックス
				my $kp = $tag[1];#名前
				for(my $i = 0; $i < 4; $i++){
					$KVEC[$ik][$i] = $tag[2+$i];#k点の整数成分(x,y,z,N)を保存
				}
				$K_NAME[$ik] = &makeEpsKname($kp);#k点の名前

				if($K_GROUP_TYPE == 1){
					$tag[6] =~ s/([CDOTS])([\w\d]*)/$1_{$2}/;
					$K_GROUP_NAME[$ik] = $tag[6];#群の名前(シェーンフリース)
				}elsif($K_GROUP_TYPE == 2){
					$K_GROUP_NAME[$ik] = $tag[7];#群の名前(国際)
				}
				$MAXDIM[$is][$ik] = 0;
			}
			next;
		}
		my ($name, $coef, $idim) = split(/\s+/, $line, 3);
		if($idim > $MAXDIM[$is][$ik]){$MAXDIM[$is][$ik] = $idim;}
		if($coef =~ /([0-9.eEdD-]+)\s*,\s*([0-9.eEdD-]+)/){
#ik,ie(エネルギー)のバンドをimr(既約表現)で簡約した係数(四捨五入した整数)、つまり、この状態に含まれる、既約表現$nameの数
			$coef = int(sqrt($1 * $1 + $2 * $2) + .5);
			if($coef > 0){
			              $nonconventional[$ik]=0;
				my ($none, $mulliken) = split(/\(/,$name,2);
				$mulliken =~ s/\)$//;
				for(my $iib = $ib; $iib<$ib+$idim; $iib++){
					if($IMR_TYPE == 0){#通常表現
						$$imr[$ik][$iib] = $none;
					}elsif($IMR_TYPE == 1){#Mulliken
						if($mulliken eq ''){
							$$imr[$ik][$iib] = $none;
							$nonconventional[$ik]=1;
						}else{
							$$imr[$ik][$iib] = $mulliken;
						}
					}
					$$imrnum[$ik][$ib] = $iimr;
				}
				$ib+=$idim;
			}
			$iimr++;
		}else{
			$$imr[$ik][$ib] = "";
			$$imrnum[$ik][$ib] = -1;
		}
   
	}
}





######################################################################
#TSPACEにおけるk点の名前をepsが表示できるギリシャ文字に変換
######################################################################
sub makeEpsKname{
    my $k = shift;
    foreach my $a (keys %TSPACE_KNAME_TABLE){#TSPACE_KNAME_TABLEの表に従ってギリシャ文字をeps用に変換
	$k =~ s/$a/$TSPACE_KNAME_TABLE{$a}/;
    }
    $k =~ s/([A-Z])P/$1'/;#二文字目のPをプライムに変化
    return $k;
}





######################################################################
#バンドの交差を表現する機能
#同じ既約表現をもつバンドをつなぎ合わせる。
#アップスピン、ダウンスピンとも機能の流れは同じ
#以下はフローチャート
#1.k-pathのループ(ikp = 0...Nkp-1、例えばikp=0がΓ-X間、ikp=2がX-L間としてk-pathの全領域をループ)
#2.ikpエリアの境界におけるikを保存(例えばikpがX-L間、X点のikが10、L点のikが20とすると、iks=10,ike=20としてk-pathの一つのエリアの境界を保存する)
#3.バンドのループ(ib = 0...NBAND-1)
#4.現在のバンドを保存(jb=ib)
#5.k=iks+1,b=ibの既約表現を保存。これがikpエリアの一つのバンドの中で既約表現が変化したかどうかを判定する基準となる。(imrcheck = imr[iks+1][ib])
#6.ikpエリア内のk点のループ(ik=iks+2...ike-1、k-pathのエリア内に存在するk点のループ(基準点ik=iks+1を除く))
#7.k=ik,b=jbの既約表現を保存。(imr2 = imr[ik][jb])
#8.imr2とimrcheckがちがうかどうか判断(違うときは9.、同じときは10.に進む)
#9.imrcheckと同じ既約表現をもち、かつk=ik-1かつb=ibのバンドと固有値が等しいバンドを探してjbを更新する(jb=lb s.t. imr[ik][lb] == imrcheck && energy[ik][lb] == energy[ik-1][jb])
#10.エネルギー配列を並べなおすため、jbを配列に保存する(ibnew[ik][ib] = jb)
#11.6.のループに戻る
#12.3.のループに戻る
#13.このk-pathエリアにおけるエネルギー配列の並べ直し(energy[ik][ib] = energy[ik][ibnew[ik][ib]], ik = iks+1...ike-1)
#14.1.のループに戻る
######################################################################
sub BandCrossSet{
	for(my $is = 0; $is < $NSPIN; $is++){
		my $imr = $IMRNAME[$is];#既約表現の番号を既約表現の比較に用いる
		my @ikpath = ();#k-pathの節のリスト
			for(my $ik = 0; $ik < $NKVEC; $ik++){
			if($LK[$is][$ik][0] == 1){
				push(@ikpath, $ik);
			}
		}
		my @pivot = ();
		for(my $ik=0;$ik<$NKVEC;$ik++){
			for(my $ib=0;$ib<$NBAND;$ib++){
			$pivot[$ik][$ib]=$ib;
			}
		}

#機能始まり、番号との対応はフローチャートを参考のこと
		for(my $ikp = 0; $ikp < (@ikpath-1); $ikp++){#1.
			my @ibnew = ();
			my ($iks, $ike) = ($ikpath[$ikp], $ikpath[$ikp+1]);#2.
#k-pathのエリア内のk点でバンドが分裂していないとこの機能は正常に作動しないので、その場合はスキップする
#			my $mdim = 0;#このk-pathエリアにおける最大縮退度
#			for(my $ik = $iks+1; $ik < $ike; $ik++){
#				if($mdim < $MAXDIM[$is][$ik]){$mdim = $MAXDIM[$is][$ik];}
#			}
#			if($mdim > 1){
#				for(my $ik = $iks; $ik<$ike+2; $ik++){
#					for(my $ib = 0;$ib <$NBAND;$ib++){
#						$pivot[$ik+1][$ib]=$pivot[$ik][$ib]
#					}
#				}
#				next ;
#			}
#この行からフローチャートに戻る
			for(my $ik = $iks+2; $ik <$ike ; $ik++){
				for(my $ib = 0; $ib <$NBAND ; $ib++){
					my $imr1 = $$imr[$ik][$ib];
					my $myimr = 0;
					my $nonflag = 1;
					for(my $iib = 0;$iib<=$ib;$iib++){#自分が何番目のその規約表現に属するバンドなのかを$myimrに入れる
						if($$imr[$ik][$iib] eq $imr1){
							$myimr ++;
						}
					}
					for(my $iib = 0; $iib < $NBAND ; $iib ++){
						my $imr3 = $$imr[$ik-1][$iib];
						if($imr1 eq $imr3){
							$myimr --;
							if($myimr == 0){
								$pivot[$ik][$ib] = $pivot[$ik-1][$iib];
								$nonflag = 0;
							}
						}
					}
					if($nonflag == 0){
						$pivot[$ik+1][$ib] = $pivot[$ik][$ib];
					}						
				}
			}
			for(my $ik = $ike-1; $ik<$ike+1; $ik++){
				for(my $ib = 0;$ib <$NBAND;$ib++){
					$pivot[$ik+1][$ib]=$pivot[$ik][$ib];
				}
			}
		}

		#
		for(my $ib = 0; $ib < $NBAND; $ib++){
			for(my $ik = 0; $ik < $NKVEC; $ik++){
				$ENERGY[$is][$ik][$pivot[$ik][$ib]]=$bufE[$is][$ik][$ib];
			}
		}

	}
}





######################################################################
#DKとLKを作成
######################################################################
sub makeDKandLK{
	$DK[0] = 0.;#最初のDKを0セット
	for(my $i=1;$i<$NKVEC;$i++) {
		my @dkc=();
		for(my $j=0;$j<3;$j++) {
			$dkc[$j] = (($KVEC[$i][$j]+0.)/$KVEC[$i][3] - ($KVEC[$i-1][$j]+0.)/$KVEC[$i-1][3]);#i+1番目のkベクトルとi番目のkベクトルの差を計算し、$dkc[0..2]に保存
		}
		$DK[$i] = $DK[$i-1] + sqrt($dkc[0]**2+$dkc[1]**2+$dkc[2]**2);#i番目のDKに、i-1番目のDKの値とdkcのノルムの和を保存する。
	}
	for(my $i=1;$i<$NKVEC;$i++){
		$DK[$i] /= $DK[$NKVEC-1];#k-pathの規格化
	}

	my @pk = (0);#kpathの節のリスト
	for(my $i = 1; $i < ($NKVEC - 1); $i++){
		if(($K_NAME[$i-1] ne $K_NAME[$i]) && ($K_NAME[$i+1] ne $K_NAME[$i])){#k点の点群が変わるところをk-pathの節と判断し、k-pathの節で既約表現を表示する
			push(@pk, $i);
		}
	}
	push(@pk, $NKVEC - 1);
	
	for(my $is = 0; $is < $NSPIN; $is++){
		for(my $ik = 0; $ik < $NKVEC; $ik++){
			my $nb = 0;#エネルギーの表示範囲内にあるバンドの数
			for(my $ib = 0; $ib < $NBAND; $ib++){
				my $en = $ENERGY[$is][$ik][$ib];
				if(($en >= $EMIN) && ($en <= $EMAX)){$nb++;}#エネルギーが表示範囲にあるとき、nbを加算
			}
			$NUM_BAND_K[$is][$ik] = $nb;
		}

		my @pkchild = ();#k-pathの節-節間のNUMIMR分割点で表示させるk点のリスト
		my @pkcenter = ();#k-pathの節-節間の中間点で表示させるk点のリスト
		for(my $ip = 1; $ip < @pk; $ip++){#pkのリストのループ
			my ($iks, $ike) = ($pk[$ip-1],$pk[$ip]);#k-pathの節-節エリアの境界におけるk点
			my $nbmax = 0;#k-pathの節-節エリアで表示されるバンドの最大数

			for(my $ik = $iks+1; $ik < $ike; $ik++){#k-pathの節-節エリアに含まれるk点のループ
				if($nbmax < $NUM_BAND_K[$is][$ik]){$nbmax = $NUM_BAND_K[$is][$ik];}
			}
			my @ikbmax = ();#表示されるバンド数がnbmaxに一致するikのリスト
			for(my $ik = $iks+1; $ik < $ike; $ik++){#k-pathの節-節エリアに含まれるk点のループ
				if($NUM_BAND_K[$is][$ik] == $nbmax){
					push(@ikbmax, $ik);
					}
			}
			my $nim=int(($DK[$ike]-$DK[$iks])/$NUMIMR);
			if($nim ==0){$nim =1;}
			for(my $num = 0;$num<$nim;$num++){
				push(@pkchild, $ikbmax[(@ikbmax+0)*($num+1)/int($nim+1)]);#ikbmax配列の分割点をpkchildに追加
			}
			push(@pkcenter, $ikbmax[(@ikbmax+0)/2]);
		}
	
#pk,pkchildに登録されたikからLKのデータをつける
		for(my $ik = 0; $ik < $NKVEC; $ik++){$LK[$is][$ik][0]=0;$LK[$is][$ik][1]=0;}
		foreach my $ik (@pk){
			$LK[$is][$ik][0] = 1;
			$LK[$is][$ik][1] = 1;
		}
		foreach my $ik (@pkchild){
			$LK[$is][$ik][0] = 2;
		}
		foreach my $ik (@pkcenter){
			$LK[$is][$ik][1] = 2;
		}
		
	}
}





######################################################################
#nfefermi.dataからフェルミ準位の値を取得する。nfefermi.dataが取得できなかった場合は0を返す。
######################################################################
sub parseFermi{
	my $arg = shift;
	my $IN;
	if(open($IN, $arg)){#nfefermi.dataを開く
		my @data = split(/\s*:\s*/,<$IN>);#nfefermi.dataの一行目からフェルミ準位(hartree)の値を取得
		close($IN);
		$arg = $data[0]*$HAR2EV;#eVに変換
	}else{
		$arg = 0.0;
	}
	return $arg;
}





######################################################################
#ヘルプメッセージの出力
######################################################################
sub help{
	print STDERR "Usage: band_symm.pl (reduce.data|reduce_up.data reduce_down.data) -erange Emin,Emax -einc dE -with_fermi nfefermi.data -linecolor (black|red|blue|...) -imrtype (default|mulliken) -imrcolor (black|red|blue|...) -font FONT,SIZE -kgrouptype (Shoenflies|HermannMauguin)\n";
}




######################################################################
#バンドクロス用subルーチン
######################################################################
