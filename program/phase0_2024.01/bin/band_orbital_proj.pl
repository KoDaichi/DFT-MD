#!/usr/bin/perl -w
#
##################################################
#
# Energy band-orbital-proj plot
# Version : 1.08
#
#  Contact address :  Phase System Consortium
#
##################################################

$eps = 1.e-4;

#################### constants
$ha2ev = 27.2113834;
$unit2ev = $ha2ev;

####################
$file_gpi="plot_band_orbproj.gnu";
$file_energyband="plot_band_orbproj.dat";

$set_erange = 'no';
$set_einc = 'no';
$set_color = 'no';
$fermi_level = 'no';
$no_dispersion = 'no';

$window_width=0.5;
$window_height=0.5;

$circle_radius=0.02;

$set_score_range = 'no';
$set_cbrange = 'no';

$set_atom_range = 'no';
$set_target_key = 'no';
$set_target_element = 'no';
$set_target_orb_l = 'no';
$set_target_orb_m = 'no';
$set_target_orb_tau = 'no';

$outfile_eps = 'orbital_projected_band.eps';
$outfile_png = 'orbital_projected_band.png';

$outfile = 'none';
$print_format = 'eps';

$vertical = 'no';

$plot_style = 1;
$threshold = 0.0;
$cbmin = 0.0;

###################
if(@ARGV<3) {
    print "[Usage] \n";
    print "    band-orbital-proj.pl EnergyDataFile KpointFile OrbProjFile -erange=Emin,Emax -einc=dE -with_fermi -atom_range=amin,amax -il=L -im=M -tau=TAU -element=X -key=I -score_range=scmin,scmax -cbrange=Cbmin,Cbmax -circle_radius=SIZE -window_width=SIZE -window_height=SIZE -vertical -color -plot_style={1,2,3} -threshold=THR -print_format={eps,png} -outfile=AAA -no_dispersion \n";
    print "\n";
    print "[Note] \n";
    print "plot_style = 1:  color is varied with orbital projection weight (default)\n";
    print "             2:  color and radius are varied \n";
    print "             3:  radius is varied \n";
    print "threshold ( lower limit ) is active when plot_style = 3 \n";
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

    } elsif($s =~/-atom_range/) {
      ($dummy,$s)=split('=',$s);
      ($atom_min,$atom_max)=split(',',$s);
      $set_atom_range = 'yes';

    } elsif($s =~/-il/) {
      ($dummy,$target_orb_l)=split('=',$s);
	$set_target_orb_l = 'yes';

    } elsif($s =~/-im/) {
      ($dummy,$target_orb_m)=split('=',$s);
	$set_target_orb_m = 'yes';

    } elsif($s =~/-tau/) {
      ($dummy,$target_orb_tau)=split('=',$s);
	$set_target_orb_tau = 'yes';

    } elsif($s =~/-element/) {
      ($dummy,$target_element)=split('=',$s);
	$set_target_element = 'yes';

    } elsif($s =~/-key/) {
      ($dummy,$target_key)=split('=',$s);
	$set_target_key = 'yes';

    } elsif($s =~/-score_range/) {
      ($dummy,$s)=split('=',$s);
      ($score_min,$score_max)=split(',',$s);
      $set_score_range = 'yes';

    } elsif($s =~/-threshold/) {
      ($dummy,$threshold)=split('=',$s);

############
    } elsif($s =~/-with_fermi/) {
        $fermi_level = 'yes';

    } elsif($s =~/-cbrange/) {
      ($dummy,$s)=split('=',$s);
      ($cbmin,$cbmax)=split(',',$s);
      $set_cbrange = 'yes';

    } elsif($s =~/-window_width/) {
      ($dummy,$window_width)=split('=',$s);
    } elsif($s =~/-window_height/) {
      ($dummy,$window_height)=split('=',$s);
    } elsif($s =~/-vertical/) {
      $vertical = 'yes';

    } elsif ($s =~/-circle_radius/) {
      ($dummy,$circle_radius)=split('=',$s);

    } elsif($s =~/-color/) {
	$set_color = 'yes';

    } elsif($s =~/-print_format/) {
      ($dummy,$print_format)=split('=',$s);
    } elsif($s =~/-outfile/) {
      ($dummy,$outfile)=split('=',$s);
    } elsif ($s =~/-plot_style/) {
      ($dummy,$plot_style)=split('=',$s);
    } elsif($s =~/-no_dispersion/) {
        $no_dispersion = 'yes';
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

###################
open(IN,$file1);
($dummy,$nkvec)=split('=',<IN>);
($dummy,$nband)=split('=',<IN>);
($dummy,$nspin)=split('=',<IN>);
($dummy,$efermi)=split('=',<IN>);
$nkvec=$nkvec/$nspin;
#print "nkvec = $nkvec\n";
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

###################
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

###################
open(IN,$file3);
<IN>;  <IN>;
($dummy,$nkvec)=split('=',<IN>);
($dummy,$nband)=split('=',<IN>);
($dummy,$nspin)=split('=',<IN>);
($dummy,$norb)=split('=',<IN>);
<IN>;<IN>;<IN>;  <IN>;

$line=<IN>;
$has_score=0;
if ( $line =~/score/ ){
    $has_score = 1;
}
#
for ( $i=0; $i <$norb; $i++ ) {
    $line = <IN>;
    @line = split(' ',$line);

    $atom = $line[1]; $il = $line[2];  $im = $line[3];  $tau = $line[4];
    $elem = $line[5]; $key=$line[6];
    if ( $has_score == 1 ) {
	$score = $line[7];
    }
    
    $Flag1 = 1;    $Flag2 = 1;	    $Flag3 = 1;	    $Flag4 = 1;
    $Flag5 = 1;    $Flag6 = 1;      $Flag7 = 1;     

    if ( $set_atom_range eq 'yes' ) {
	if ( $atom < $atom_min || $atom > $atom_max ) {
	    $Flag1 = 0;
	}
    }
    if ( $set_target_orb_l eq 'yes' && $il != $target_orb_l ) {
	$Flag2 = 0;
    }
    if ( $set_target_orb_m eq 'yes' && $im != $target_orb_m ) {
	$Flag3 = 0;
    }
    if ( $set_target_orb_tau eq 'yes' && $tau != $target_orb_tau ) {
	$Flag4 = 0;
    }
    if ( $set_target_element eq 'yes' && $elem ne $target_element ) {
	$Flag5 = 0;
    }
    if ( $set_target_key eq 'yes' && $key ne $target_key ) {
	$Flag6 = 0;
    }
    if ( $has_score == 1 ) {
	if ( $set_score_range eq 'yes' ) {
	    if ( $score < $score_min || $score > $score_max ) {
		$Flag7 = 0;
	    }
	}
    }

    if ( $Flag1 *$Flag2 *$Flag3 *$Flag4 *$Flag5 *$Flag6 *$Flag7 == 1 ) {
	$sum_flag[$i] = 1;
    } else {
	$sum_flag[$i] = 0;
    }	
}

<IN>;  <IN>;

$nkvec=$nkvec/$nspin;

$i=0;
while($i < $nkvec) {
    for ( $s=0; $s<$nspin; $s++ ) {
	for ( $nn=0; $nn <$nband; $nn++ ) {
	    ${$weight[$s][$i]}[$nn] = 0.0;
	}
    }
    $i++;
}

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
#
	    for ( $ll=0; $ll <$norb; $ll++ ){
		$line = <IN>;
		if ( $line =~/ia/ ) {
#		    @line3 = split(' ',$line);
#    		    $iorb = $line3[0];
		    
		    @w=();
		    while(@w < $nband) {
			push @w, split(' ',<IN>);
		    }

#		    if ( $sum_flag[$iorb-1] eq 1 ) {
		    if ( $sum_flag[$ll] eq 1 ) {
			for ( $nn=0; $nn <$nband; $nn++ ) {
			    ${$weight[$s][$i]}[$nn] = ${$weight[$s][$i]}[$nn] +$w[$nn];
			}
		    }
		}
	    }
	    <IN>;    <IN>;
	}
    } 
    $i++;
}
close(IN);


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

#print OUT "set size $window_width,1.0\n";

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
print OUT "set gr \n";

print OUT "emin=$emin \n";
print OUT "emax=$emax \n";

print OUT "set xr [0.00:$dk[$nkvec-1]] \n";
print OUT "set yr [emin\:emax] \n";

if($set_cbrange eq 'yes') {
    print OUT "cbmin=$cbmin \n";
    print OUT "cbmax=$cbmax \n";
    print OUT "set cbrange [cbmin:cbmax] \n";
}
if ($set_color eq 'no') {
    print OUT "set palette defined ( 0 'white', 1 'gray30' ) \n";
    print OUT "set style line 100 lc palette \n";
    print OUT "set style line 200 lc rgb \"black\" \n";
} else {
    print OUT "set palette rgb 22,13,-31 \n";
    print OUT "set style line 100 lc palette \n";
    print OUT "set style line 200 lc rgb \"web-blue\" \n";
}

if ($fermi_level eq 'yes') {
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

if ( $nspin == 1 ) {
    if ( $plot_style == 3 ) {print OUT "threshold=$threshold \n";}

    print OUT "circle_radius=$circle_radius \n";

    if ( $plot_style == 1 ){
	print OUT "plot \"$file_energyband\" u 1:2:(circle_radius):(\$3) w circles fill transparent solid 0.2 noborder ls 100";

    } elsif ( $plot_style == 2 ) {
	print OUT "plot \"$file_energyband\" u 1:2:(circle_radius*\$3):(\$3) w circles fill transparent solid 0.2 noborder ls 100";

    } elsif ( $plot_style == 3 ) {
	print OUT "threshold=$threshold \n";
	print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.2 noborder ls 200";
    }
    if ( $no_dispersion eq 'yes' ){
	print OUT "\n";
    } else {
	print OUT ", \\\n";
	print OUT "\"$file_energyband\" u 1:2 w l lw 0.2 lc \"black\" \n";
    }

} else {
    if ( $plot_style == 3 ) {print OUT "threshold=$threshold \n";}

    print OUT "circle_radius=$circle_radius \n";
    print OUT "set multiplot layout $lay_x, $lay_y \n";

    if ( $plot_style == 1 ){
	print OUT "plot \"$file_energyband\" u 1:2:(circle_radius):(\$3) w circles fill transparent solid 0.8 noborder ls 100";
    } elsif ( $plot_style == 2 ) {
	print OUT "plot \"$file_energyband\" u 1:2:(circle_radius*\$3):(\$3) w circles fill transparent solid 0.8 noborder ls 100";
    } elsif ( $plot_style == 3 ) {
	print OUT "plot \"$file_energyband\" u 1:(\$3>=threshold ? \$2:1/0):(circle_radius*(\$3-threshold)) w circles fill transparent solid 0.8 noborder ls 200";
    }

    if ( $no_dispersion eq 'yes' ){
	print OUT "\n";
    } else {
	print OUT ", \\\n";
	print OUT "\"$file_energyband\" u 1:2 w l lw 0.2 lc \"black\" \n";
    }

    if ( $plot_style == 1 ){
	print OUT "plot \"$file_energyband\" u 1:4:(circle_radius):(\$5) w circles fill transparent solid 0.8 noborder ls 100";
    } elsif ( $plot_style == 2 ) {
	print OUT "plot \"$file_energyband\" u 1:4:(circle_radius*\$5):(\$5) w circles fill transparent solid 0.8 noborder ls 100";
    } elsif ( $plot_style == 3 ) {
	print OUT "plot \"$file_energyband\" u 1:(\$5>=threshold ? \$4:1/0):(circle_radius*(\$5-threshold)) w circles fill transparent solid 0.8 noborder ls 200";
    }

    if ( $no_dispersion eq 'yes' ){
	print OUT "\n";
    } else {
	print OUT ", \\\n";
	print OUT "\"$file_energyband\" u 1:4 w l lw 0.2 lc \"black\" \n";
    }

    print OUT "unset multiplot \n";
}

close(OUT);

######################################
$nstate = @{$energy[0][0]};

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

if ( $print_format eq "eps" or $print_format eq "png" ) {
    `gnuplot $file_gpi`;
}
