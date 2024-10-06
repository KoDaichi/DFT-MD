#!/usr/bin/perl -w
#
# Energy band-unfold plot
# Version : 1.07
#
#  by ASMS
#
#  Contact address :  Phase System Consortium
#

$eps = 1.e-4;

#################### constants
$ha2ev = 27.2113834;
$unit2ev = $ha2ev;

$PAI=3.1415926535e0;

####################
$file_gpi="plot_band_unfolding.gnu";
$file_energyband="plot_band_energy.dat";
$file_energymap="plot_band_energy.map";

$set_erange = 'no';
$set_einc = 'no';
$window_width=0.5;
$window_height=0.5;

#### for plotting spectral weight
$circle_radius=0.02;
$with_dispersion = 'no';
$fermi_level = 'no';

#### for plotting spectral function
$plot_spectral_func = 'no';
$set_color = 'no';
$set_cbrange = 'no';
$ndiv = 400;
$sigma = 0.05e0;        
$line_width = 4;

$outfile_eps = 'unfolded_band.eps';
$outfile_png = 'unfolded_band.png';

$outfile = 'none';
$print_format = 'eps';

$vertical = 'no';

$threshold = 0.0;

###################
if(@ARGV<3) {
    print "[Usage]\n";
    print "* Plotting spectral weight \n";
    print "    band-unfold.pl EnergyDataFile KpointFile SpectralWeightFile -erange=Emin,Emax -einc=dE -with_fermi -threshold=THR -window_width=SIZE -window_height=SIZE -vertical -with_dispersion -circle_radius=SIZE -color -print_format={eps,png} -outfile=AAA \n";
    print "\n";
    print "* Plotting spectral function \n";
    print "    band-unfold.pl EnergyDataFile KpointFile SpectralWeightFile -spectral_func -erange=Emin,Emax -einc=dE -with_fermi -window_width=SIZE -window_height=SIZE -vertical -ndiv=VAL -sigma=VAL -line_width=VAL -cbrange=Cbmin,Cbmax -color -print_format={eps,png} -outfile=AAA \n";
    print "\n";
  exit;

} else {
  foreach $s ( @ARGV ) {
    if($s =~/-erange/) {
      ($dummy,$s)=split('=',$s);
      ($emin,$emax)=split(',',$s);
      $set_erange = 'yes';
    } elsif ($s =~/-einc/) {
      ($dummy,$einc)=split('=',$s);
	$set_einc = 'yes';
    } elsif($s =~/-window_width/) {
      ($dummy,$window_width)=split('=',$s);
    } elsif($s =~/-window_height/) {
      ($dummy,$window_height)=split('=',$s);
    } elsif($s =~/-vertical/) {
      $vertical = 'yes';

#### for plotting spectral weight
    } elsif($s =~/-circle_radius/) {
      ($dummy,$circle_radius)=split('=',$s);
    } elsif($s =~/-with_fermi/) {
	$fermi_level = 'yes';
    } elsif($s =~/-threshold/) {
      ($dummy,$threshold)=split('=',$s);
    } elsif($s =~/-with_dispersion/) {
	$with_dispersion = 'yes';

#### for plotting spectral function
    } elsif($s =~/-spectral_func/) {
      $plot_spectral_func = 'yes';
    } elsif($s =~/-ndiv/) {
      ($dummy,$ndiv)=split('=',$s);
    } elsif($s =~/-sigma/) {
      ($dummy,$sigma)=split('=',$s);
    } elsif($s =~/-line_width/) {
      ($dummy,$line_width)=split('=',$s);
    } elsif($s =~/-cbrange/) {
      ($dummy,$s)=split('=',$s);
      ($cbmin,$cbmax)=split(',',$s);
      $set_cbrange = 'yes';
    } elsif($s =~/-line_width/) {
      ($dummy,$line_width)=split('=',$s);
    } elsif($s =~/-color/) {
      $set_color = 'yes';
#####
      
    } elsif($s =~/-print_format/) {
      ($dummy,$print_format)=split('=',$s);
    } elsif($s =~/-outfile/) {
      ($dummy,$outfile)=split('=',$s);
    }
  }
}
if ( $print_format eq 'eps' ) {
    if ( $outfile eq 'none' ) {
        $outfile = $outfile_eps;
    }
}
if ( $print_format eq 'png' ) {
    if ( $outfile eq 'none' ) {
        $outfile = $outfile_png;
    }
}

$file1 = $ARGV[0];
$file2 = $ARGV[1];
$file3 = $ARGV[2];

######################################
open(IN,$file1);
($dummy,$nkvec)=split('=',<IN>);
($dummy,$nband)=split('=',<IN>);
($dummy,$nspin)=split('=',<IN>);
($dummy,$efermi)=split('=',<IN>);
$nkvec=$nkvec/$nspin;

$i=0;
while($i < $nkvec) {
    $line=<IN>;
    if($line=~/ ===/) {
	for($s=0;$s<$nspin;$s++) {
	    @line = split('\(',<IN>);
	    @line = split('\)',$line[1]);
	    @line = split(' ',$line[0]);
	    $kvec[$i][0]=$line[0];
	    $kvec[$i][1]=$line[1];
	    $kvec[$i][2]=$line[2];

	    @e=();
	    while(@e < $nband) {
		push @e, split(' ',<IN>);
	    }
	    @{$energy[$s][$i]} = @e; 
	    if($nspin == 2 && $s == 0) {
		<IN>;
	    }
	}
	$i++;
    } 
}
close(IN);

######################################
open(IN,$file2);
<IN>;
@{$rlvec[0]} = split(' ',<IN>);
@{$rlvec[1]} = split(' ',<IN>);
@{$rlvec[2]} = split(' ',<IN>);
$i=0;
while($line=<IN>) {
  ($line,$spklabel[$i]) = split('#',$line);
  chomp($spklabel[$i]);
  ($n1,$n2,$n3,$nd) = split(' ',$line);
  $spk[$i][0] = $n1/$nd;
  $spk[$i][1] = $n2/$nd;
  $spk[$i][2] = $n3/$nd;
  $i++;
}
close(IN);

######################################
open(IN,$file3);
($dummy,$nkvec)=split('=',<IN>);
($dummy,$nband)=split('=',<IN>);
($dummy,$nspin)=split('=',<IN>);
$line=<IN>;
$nkvec=$nkvec/$nspin;

$i=0;
while($i < $nkvec) {
    for($s=0;$s<$nspin;$s++) {
	$line=<IN>;
	if($line=~/ik =/) {
	    @line2 = split('\(',$line);
	    @line2 = split('\)',$line2[1]);
	    @line2 = split(' ',$line2[0]);
	    $kvec[$i][0]=$line2[0];
	    $kvec[$i][1]=$line2[1];
	    $kvec[$i][2]=$line2[2];
	    
	    @w=();
	    while(@w < $nband) {
		push @w, split(' ',<IN>);
	    }
	    @{$weight[$s][$i]} = @w; 
	    $line = <IN>;
	}
    } 
    $i++;
}
close(IN);

###############################
for($k=0;$k<$nkvec;$k++) {
  $label[$k] = 'none';
  for($i=0;$i<@spk;$i++) {
    if(abs($spk[$i][0]-$kvec[$k][0]) < $eps && 
       abs($spk[$i][1]-$kvec[$k][1]) < $eps &&
       abs($spk[$i][2]-$kvec[$k][2]) < $eps ) {
      $label[$k] = $spklabel[$i];
      last;
    }
  }
  for($s=0;$s<$nspin;$s++) {
    for($i=0;$i<@{$energy[$s][$k]};$i++) {
      $energy[$s][$k][$i] -= $efermi;
      $energy[$s][$k][$i] *= $unit2ev;
    }
  }
}

if($set_erange eq 'no') {
  $emin=1e10;
  $emax=-1e10;
  for($s=0;$s<$nspin;$s++) {
    for($k=0;$k<$nkvec;$k++) {
      for($i=0;$i<@{$energy[$s][$k]};$i++) {
        if($energy[$s][$k][$i] < $emin) {
          $emin = $energy[$s][$k][$i];
	  } elsif($energy[$s][$k][$i] > $emax) {
	    $emax = $energy[$s][$k][$i];
	  }
	}
    }
  }
  $shift=0;
  $shift=-5 if($emin<0);
  $emin = int($emin/5)*5 + $shift;
  $shift= 5 if($emax>0);
  $emax = int($emax/5)*5 + $shift;
}

@dk = ( 0.0, 0.0, 0.0 );
for($i=1;$i<$nkvec;$i++) {
  @dkc=();
  for($j=0;$j<3;$j++) {
    for($k=0;$k<3;$k++) {
      $dkc[$j] += $rlvec[$j][$k]*($kvec[$i][$k]-$kvec[$i-1][$k]);
    }
  }
  $dk[$i] = $dk[$i-1] + sqrt($dkc[0]**2+$dkc[1]**2+$dkc[2]**2);
}

##################################
open(OUT,">$file_gpi");

$wscale_x = 2.0 *$window_width;
$wscale_y = 2.0 *$window_height;

if ( $nspin == 1 ) {
    $wscale_x = $wscale_x *0.5;
}

if ( $print_format eq 'eps' ) {
    $w1 = 5;  $w2 = 3.5;
    print OUT "wsx=$wscale_x \n";
    print OUT "wsy=$wscale_y \n";
    if($set_color eq 'no') {
        print OUT "set terminal postscript eps enhanced solid";
    } else {
        print OUT "set terminal postscript eps enhanced color solid";
    }
    print OUT " size $w1*wsx, $w2*wsy \n";

} elsif ( $print_format eq 'png' ) {
    $w1 = 640;  $w2 = 480;
    print OUT "wsx=$wscale_x \n";
    print OUT "wsy=$wscale_y \n";
    print OUT "set terminal pngcairo enhanced ";
    print OUT " size $w1*wsx, $w2*wsy \n";
}

print OUT "set output \"$outfile\" \n";

print OUT "set nokey\n";
if ( $plot_spectral_func eq 'yes' )  {
    print OUT "set view 0 \, 0 \n";
    print OUT "set pm3d map \n";
}
print OUT "emin=$emin \n";
print OUT "emax=$emax \n";
 
print OUT "set xr [0.00:$dk[$nkvec-1]] \n";
print OUT "set yr [emin\:emax] \n";

if ($set_cbrange eq 'yes') {
    print OUT "cbmin=$cbmin \n";
    print OUT "cbmax=$cbmax \n";
    print OUT "set cbrange [cbmin\:cbmax] \n";
}
if ( $plot_spectral_func eq 'yes' )  {
    if ($set_color eq 'no') {
	print OUT "set palette defined ( 0 'white', 1 'light-gray', 2 'dark-gray', 3 'gray50', 4 'gray30' ) \n";
    } else {
	print OUT "set palette defined ( 0 'white', 1 'light-blue', 2 'skyblue', 3 'royalblue', 4 'medium-blue', 5 'midnight-blue' ) \n";
    }
} else {
    if ($set_color eq 'no') {
	print OUT "set style line 100 lc rgb 'gray30' \n";
    } else {
	print OUT "set style line 100 lc rgb 'red' \n";
#	print OUT "set style line 100 lc rgb 'web-blue' \n";
    }
}

if($fermi_level eq 'yes') {
    print OUT "set arrow front from 0.0, 0,0 to $dk[$nkvec-1], 0.0 nohead lt 0.2 \n";
}
for($i=1;$i<@label-1;$i++) {
  if(! (lc($label[$i]) =~ 'none')) {
    print OUT "set arrow front from $dk[$i],emin to $dk[$i],emax nohead lt -1\n";
  }
}

print OUT "set ylabel \"Energy (eV)\"\n";
print OUT "set xtics ( ";
for($i=0;$i<@label;$i++) {
  if(! (lc($label[$i]) =~ 'none')) {
    print OUT "\"$label[$i]\" $dk[$i]";
    print OUT ", "if($i != $nkvec-1);
  }
}
print OUT ")\n";

if($set_einc eq 'yes') {
  print OUT "set ytics $einc\n";
}

if ( $vertical eq 'yes' ) {
    $lay_x = 2;  $lay_y = 1;
} else {
    $lay_x = 1;  $lay_y = 2;
}

if ( $plot_spectral_func eq 'yes' )  {
    if ( $nspin == 1 ) {
	if ( $with_dispersion eq 'yes' ) {
	    print OUT "plot \"$file_energymap\" u 1:2:3 w image, ";
	    print OUT "\"$file_energyband\" u 1:2 w l lw 0.1 lc \"black\" \n";
	} else {
#	    print OUT "splot \"$file_energymap\" u 1:2:3 w l lt 1 lw $line_width lc palette \n";
	    print OUT "plot \"$file_energymap\" u 1:2:3 w image \n";
	}
    } else {
	if ( $with_dispersion eq 'yes' ) {
	    print OUT "set multiplot layout $lay_x, $lay_y \n";
	    print OUT "plot \"$file_energymap\" u 1:2:3 w image, ";
	    print OUT "\"$file_energyband\" u 1:2 w l lw 0.1 lc \"black\" \n";
	    print OUT "plot \"$file_energymap\" u 1:2:4 w image, ";
	    print OUT "\"$file_energyband\" u 1:4 w l lw 0.1 lc \"black\" \n";
	    print OUT "unset multiplot \n";
	} else {
	    print OUT "set multiplot layout $lay_x, $lay_y \n";
#	    print OUT "splot \"$file_energymap\" u 1:2:3 w l lt 1 lw $line_width lc palette \n";
#	    print OUT "splot \"$file_energymap\" u 1:2:4 w l lt 2 lw $line_width lc palette \n";
	    print OUT "plot \"$file_energymap\" u 1:2:3 w image \n";
	    print OUT "plot \"$file_energymap\" u 1:2:4 w image \n";
	    print OUT "unset multiplot \n";
	}
    }
} else {
    print OUT "circle_radius=$circle_radius \n";
    print OUT "threshold=$threshold \n";
    if ( $nspin == 1 ) {
	if ( $with_dispersion eq 'yes' ) {
	    print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.8 noborder ls 100, ";
	    print OUT "\"$file_energyband\" u 1:2 w l lw 0.1 lc \"black\" \n";
	} else {
	    print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.8 noborder ls 100 \n";
	}
    } else {
	print OUT "set multiplot layout $lay_x, $lay_y \n";
	if ( $with_dispersion eq 'yes' ) {
	    print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.8 noborder ls 100, ";
	    print OUT "\"$file_energyband\" u 1:2 w l lw 0.1 lc \"black\" \n";
	    print OUT "plot \"$file_energyband\" u 1:(\$5>=threshold ? \$4:1/0):(circle_radius*(\$5-threshold)) w circles fill transparent solid 0.8 noborder ls 100, ";
	    print OUT "\"$file_energyband\" u 1:4 w l lw 0.1 lc \"black\" \n";
	} else {
	    print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.8 noborder ls 100 \n ";
	    print OUT "plot \"$file_energyband\" u 1:(\$5>=threshold ? \$4:1/0):(circle_radius*(\$5-threshold)) w circles fill transparent solid 0.8 noborder ls 100 \n ";
	}
	print OUT "unset multiplot \n";
    }
}
close(OUT);

######################################o
$nstate = @{$energy[0][0]};

if ( $plot_spectral_func  eq 'yes' )  {
    open(OUT,">$file_energymap");

    $ediv = ( $emax -$emin )/ $ndiv;
    $coeff = 1.0 /sqrt( 2.0 *$PAI ) /$sigma *$ediv;

    for($i=0;$i<$nkvec;$i++) {
	for ( $nn=0; $nn <$ndiv; $nn++ ) {
	    $ene = $emin +$ediv *$nn;
	    for ($s=0; $s<$nspin; $s++ ) {
		$csum[$s] = 0.0;
	    }
	    for ($s=0; $s<$nspin; $s++ ) {
		for ( $n=0; $n <$nstate; $n++ ) {
		    $ee = ${$energy[$s][$i]}[$n];
		    $ww = ${$weight[$s][$i]}[$n];
		    $csum[$s] = $csum[$s] +$ww *exp( -( $ene-$ee )**2 /2.0 /$sigma**2 );
		}
		$csum[$s] = $csum[$s] *$coeff;
	    }
	    if ( $nspin == 1 ) {
		print OUT "$dk[$i] $ene $csum[0] \n";
	    } else {
		print OUT "$dk[$i] $ene $csum[0] $csum[1] \n";
	    }
	}
	print OUT "\n";
    }
    print OUT "end\n";
    close(OUT);
}
######################################o
if ( $plot_spectral_func eq 'no' )  {
    open(OUT,">$file_energyband");

    for ( $n=0; $n <$nstate; $n++ ) {
	for($i=0;$i<$nkvec;$i++) {
	    if ( $nspin == 1 ) {
		$ee1 = ${$energy[0][$i]}[$n];
		$ww1 = ${$weight[0][$i]}[$n];
		print OUT "$dk[$i] $ee1 $ww1 \n";
	    } else {
		$ee1 = ${$energy[0][$i]}[$n];
		$ww1 = ${$weight[0][$i]}[$n];
		$ee2 = ${$energy[1][$i]}[$n];
		$ww2 = ${$weight[1][$i]}[$n];
		print OUT "$dk[$i] $ee1 $ww1 $ee2 $ww2 \n";
	    }
	}
	print OUT "\n";
    }
    print OUT "end\n";
    close(OUT);
}

if ( $print_format eq "eps" or $print_format eq "png" ) {
    `gnuplot $file_gpi`;
}
