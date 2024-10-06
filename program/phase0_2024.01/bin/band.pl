#!/usr/bin/perl -w
#
# Energy band plot
# Version : 1.04
#
# Copyright (C) 2003-2004, Takenori YAMAMOTO
# Institute of Industrial Science, The University of Tokyo
#

$eps = 1.e-4;
$psize = 0.5;

# constants
$ha2ev = 27.2113834;
$unit2ev = $ha2ev;

$file_energyband = "energyband_$$.dat";
$file_gpi = "bandplot_$$.gpi";

$plot_type = 'lines'; # {solid_circles|lines}
$set_erange = 'no';
$fermi_level = 'no';
$set_einc = 'no';
$set_color = 'no';

$width=0.5;

if(@ARGV<2) {
  print "Usage: band.pl EnergyDataFile KpointFile -erange=Emin,Emax -einc=dE -ptype={solid_circles|lines} -with_fermi -width=SIZE -color\n";
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
    } elsif($s =~/-ptype/) {
      ($dummy,$plot_type)=split('=',$s);
    } elsif($s =~/-with_fermi/) {
      $fermi_level = 'yes';
    } elsif($s =~/-width/) {
      ($dummy,$width)=split('=',$s);
    } elsif($s =~/-color/) {
      $set_color = 'yes';
    }
  }
}

$file1 = $ARGV[0];
$file2 = $ARGV[1];

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
#        <IN>;<IN>;
        <IN>;
	}
      #print "$s $i : @{$kvec[$i]} : @e\n";
    }
    $i++;
  } 
}
close(IN);

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
  #print "$label[$k] @{$kvec[$k]}\n";
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

open(OUT,">$file_gpi");
if($set_color eq 'no') {
  print OUT "set terminal postscript eps enhanced solid\n";
} else {
  print OUT "set terminal postscript eps enhanced color solid\n";
}
print OUT "set size $width,1.0\n";
print OUT "set output \"band_structure.eps\"\n";
print OUT "set nokey\n";
print OUT "set ylabel \"Energy (eV)\"\n";

for($i=1;$i<@label-1;$i++) {
  if(! (lc($label[$i]) =~ 'none')) {
    print OUT "set arrow from $dk[$i],$emin to $dk[$i],$emax nohead lt -1\n";
  }
}

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

$nstate = @{$energy[0][0]};
if($plot_type =~/^p/ or $plot_type =~/^s/) {
  if($plot_type =~/^p/) {
    $point_type = 0;
    $point_type_down = 0;
  } else {
    if($nspin == 1) {
      $point_type = 7;
    } else {
      $point_type = 7;
      $point_type_down = 1;
      $point_type_down = 7 if($set_color eq 'yes');
    }
  }
  print OUT "plot [][$emin\:$emax] \\\n";
  if($fermi_level eq 'yes') {
    print OUT "0 w l , \\\n";
  }
  for($i=0;$i<$nstate-1;$i++) {
    $n=$i+2;
    print OUT "\"$file_energyband\" u 1:$n w p lt 1 pt $point_type ps $psize, \\\n";
  }
  $n=$nstate+1;
  print OUT "\"$file_energyband\" u 1:$n w p lt 1 pt $point_type ps $psize";
  if($nspin > 1) {
    print OUT ", \\\n";
  } else {
    print OUT "\n";
  }
  if($nspin > 1) {
    if($set_color eq 'yes') {
      $linetype = 3;
    } else {
      $linetype = 1;
    }
    for($i=0;$i<$nstate-1;$i++) {
      $n=$nstate+$i+2;
      print OUT "\"$file_energyband\" u 1:$n w p lt $linetype pt $point_type_down ps $psize, \\\n";
    }
    $n=2*$nstate+1;
    print OUT "\"$file_energyband\" u 1:$n w p lt $linetype pt $point_type_down ps $psize\n";
  }
} elsif($plot_type =~/^l/) {
  print OUT "plot [][$emin\:$emax] \\\n";
  if($fermi_level eq 'yes') {
    print OUT "0 w l lt 0, \\\n";
  }
  for($i=0;$i<$nstate-1;$i++) {
    $n=$i+2;
    print OUT "\"$file_energyband\" u 1:$n w l lt 1 lw $psize, \\\n";
  }
  $n=$nstate+1;
  print OUT "\"$file_energyband\" u 1:$n w l lt 1 lw $psize";
  if($nspin == 1) {
    print OUT "\n";
  } else {
    print OUT ", \\\n";
  }
  if($nspin > 1) {
    if($set_color eq 'yes') {
      $linetype = 3;
    } else {
      $linetype = 0;
    }
    for($i=0;$i<$nstate-1;$i++) {
      $n=$nstate+$i+2;
      print OUT "\"$file_energyband\" u 1:$n w l lt $linetype lw $psize, \\\n";
    }
    $n=2*$nstate+1;
    print OUT "\"$file_energyband\" u 1:$n w l lt $linetype lw $psize\n";
  }

} else {
  #print "Set -ptype={points|solid_circles|curves}.\n";
  print "Set -ptype={solid_circles|lines}.\n";
  exit 1;
}
close(OUT);

open(OUT,">$file_energyband");
for($i=0;$i<$nkvec;$i++) {
  if($nspin == 1) {
    print OUT "$dk[$i] @{$energy[0][$i]}\n";
  } else {
    print OUT "$dk[$i] @{$energy[0][$i]} @{$energy[1][$i]}\n";
  }
}
print OUT "end\n";
close(OUT);

`gnuplot $file_gpi`;
`rm -f $file_energyband`;
`rm -f $file_gpi`;
