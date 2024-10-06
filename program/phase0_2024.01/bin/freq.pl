#!/usr/bin/perl
#
# A visualizaion program for vibrational freqencies
# Version : 1.01
#
# Copyright (C) 2004, Takenori YAMAMOTO
# Institute of Industrial Science, The University of Tokyo
#

if(@ARGV == 0) {
   print "*** A visualization program for vibrational freqencies ***\n";
   print "Usage: freq.pl [-width=W] [-height=H] [-nrep=N] {-solid|-mol|-ignored_modes=LIST} mode.data\n";
   exit;
}

@ignored_modes = ();
$ignored_trans = 'yes';
$ignored_rot = 'no';
$set_nrep = 'no';
$width = 1;
$height = 1;
while( @ARGV > 1) {
  $opt = shift @ARGV;
  if($opt eq '-solid') {
     $ignored_trans = 'yes';
  } elsif($opt eq '-mol') {
     $ignored_trans = 'yes';
     $ignored_rot = 'yes';
  } elsif($opt =~/-ignored_modes/) {
     ($dummy,$list)=split('=',$opt);
     @ignored_modes=split(' ',$list);
     ##$ignored_trans = 'no';
     $ignored_rot = 'no';
  } elsif($opt =~/-nrep/) {
     ($dummy,$nrep)=split('=',$opt);
     $set_nrep = 'yes';
  } elsif($opt =~/-width/) {
     ($dummy,$width)=split('=',$opt);
  } elsif($opt =~/-height/) {
     ($dummy,$height)=split('=',$opt);
  }
}
$file_modes = shift @ARGV;

open(IN,"$file_modes");
$dielectric_mode = 'no';
while($line=<IN>) {
  if($line=~/Mode effective charge/) {
    $dielectric_mode = 'yes';
  }
}
close(IN);

open(IN,"$file_modes");
<IN>;
@{$a[0]} = split(' ',<IN>);
@{$a[1]} = split(' ',<IN>);
@{$a[2]} = split(' ',<IN>);
<IN>;
($dummy,$natom) = split('=',<IN>);
for($i=0;$i<$natom;$i++) {
   ($no,$pos[$i][0],$pos[$i][1],$pos[$i][2],$mass[$i],$name[$i]) = split(' ',<IN>);
   print "$no @{$pos[$i]} $mass[$i]\n";
}
<IN>;
($dummy,$nmode,$dummy,$natom) = split(' ',<IN>);
$max_norm = 0;
for($j=0;$j<$nmode;$j++) {
  print "mode: $j\n";
  ($dummy,$dummy,$rep[$j],$active[$j])=split(' ',<IN>);
  ($dummy,$omega[$j],$dummy,$dummy,$omega_ev[$j],$dummy,$dummy,$freq[$j]) = split(' ',<IN>);
  $omega[$j] =~ s/D/E/;
  $omega_ev[$j] =~ s/D/E/;
  $freq[$j] =~ s/D/E/;
  for($i=0;$i<$natom;$i++) {
    ($no,@{$u[$j][$i]}) = split(' ',<IN>);
    $norm[$j][$i] = 0;
    for($k=0;$k<3;$k++) {
      $im = $i*3+$k;
      $xi[$j][$im] = $u[$j][$i][$k];
      $u[$j][$i][$k] = $u[$j][$i][$k]/sqrt($mass[$i]);
      $norm[$j][$i] += $u[$j][$i][$k]**2;
    }
    $norm[$j][$i] = sqrt($norm[$j][$i]);
    if($norm[$j][$i] > $max_norm) {
      $max_norm = $norm[$j][$i];
    }
    print "$no  u= @{$u[$j][$i]} norm= $norm[$j][$i]\n";
  }
  if($dielectric_mode eq 'yes') {
   <IN>; # label
   <IN>; # mode effective charge
  }
}

# set ignored modes
if($ignored_trans eq 'yes') {
  for($i=0;$i<3;$i++) {
    for($id=0;$id<$nmode;$id++) {
      $xi_trans[$id] = 0.0;
    }
    $norm = 0.0;
    for($j=0;$j<$natom;$j++) {
      $id = $j*3+$i;
      $xi_trans[$id] = sqrt($mass[$j]);
      $norm += $xi_trans[$id]**2;
    }
    $norm = 1.0/sqrt($norm);
    for($id=0;$id<$nmode;$id++) {
      $xi_trans[$id] = $xi_trans[$id]*$norm; 
    }
    for($im=0;$im<$nmode;$im++) {
      $ctrans = 0.0;
      for($id=0;$id<$nmode;$id++) { 
        $ctrans += $xi_trans[$id]*$xi[$im][$id];
      }
      #print "i=$i im=$im ctrans=$ctrans\n";
      if(abs($ctrans) > 0.01) {
        push @ignored_modes, $im;
      }
    }
  }
}
if($ignored_rot eq 'yes') {
  for($i=0;$i<3;$i++) {
    for($id=0;$id<$nmode;$id++) {
      $xi_rot[$id] = 0.0;
    }
    $norm = 0.0;
    for($j=0;$j<$natom;$j++) {
      $i1 = $i%3;
      $i2 = ($i+1)%3;
      $id1 = $j*3+$i1;
      $id2 = $j*3+$i2;
      $xi_rot[$id2] = sqrt($mass[$j])*$pos[$j][$i1];
      $xi_rot[$id1] = -sqrt($mass[$j])*$pos[$j][$i2];
      $norm += $xi_rot[$id1]**2 + $xi_rot[$id2]**2;
    }
    $norm = 1.0/sqrt($norm);
    for($id=0;$id<$nmode;$id++) {
      $xi_rot[$id] = $xi_rot[$id]*$norm; 
    }
    for($im=0;$im<$nmode;$im++) {
      $crot = 0.0;
      for($id=0;$id<$nmode;$id++) { 
        $crot += $xi_rot[$id]*$xi[$im][$id];
      }
      #print "i=$i im=$im crot=$crot\n";
      if(abs($crot) > 0.01) {
        push @ignored_modes, $im;
      }
    }
  }
}
print "Ignored modes: @ignored_modes\n";

# Representations
$im0=0;
L1: for($im=0;$im<$nmode;$im++) {
  $yn = 'no';
  for($ig=0;$ig<@ignored_modes;$ig++) {
    if($ignored_modes[$ig] == $im) {
      $yn = 'yes';
    }
  }
  if($yn eq 'no') {
    $im0 = $im;
    last L1;
  }
}
#print "im0 = $im0\n";

@rep_names = ( $rep[$im0] );
@act_names = ( $active[$im0] );
L2: for($j=$im0+1;$j<$nmode;$j++) {
  for($ig=0;$ig<@ignored_modes;$ig++) {
    next L2 if($ignored_modes[$ig] == $j);
  }
  $match = 'no';
  for($i=0;$i<@rep_names;$i++) {
    if($rep_names[$i] eq $rep[$j]) {
      $match = 'yes';
    }
  }
  if($match eq 'no') {
    push @rep_names, $rep[$j];
    push @act_names, $active[$j];
  }
}
print "@rep_names\n";
print "@act_names\n";

# Frequencies
for($i=0;$i<@rep_names;$i++) {
  $k = 0;
  L3: for($j=$im0;$j<$nmode;$j++) {
    for($ig=0;$ig<@ignored_modes;$ig++) {
      next L3 if($ignored_modes[$ig] == $j);
    }
    if($rep_names[$i] eq $rep[$j]) {
      $freq_rep[$i][$k] = $freq[$j];
      $k++;
    }
  }
  print "$rep_names[$i] @{$freq_rep[$i]}\n";
  $ndeg = 1;
  if($rep_names[$i] =~/^E/) {
    $ndeg = 2; 
  } elsif($rep_names[$i] =~/^T/) {
    $ndeg = 3; 
  }
  if($ndeg > 1) {
    @freq_tmp = ();
    push @freq_tmp, $freq_rep[$i][0];
    for($j=1;$j<@{$freq_rep[$i]};$j++) {
      $match = 'no';
      for($k=0;$k<@freq_tmp;$k++) {
        if(abs($freq_rep[$i][$j]-$freq_tmp[$k]) < 0.1) {
          $match = 'yes';
        }
      }
      if($match eq 'no') {
        push @freq_tmp, $freq_rep[$i][$j];
      }
    }
    @{$freq_rep[$i]} = ();
    for($j=0;$j<@freq_tmp;$j++) {
      $freq_rep[$i][$j] = $freq_tmp[$j];
    }
    print " ==> $rep_names[$i] @{$freq_rep[$i]}\n";
  }
}

$min_freq = $freq[$im0];
$max_freq = $freq[$im0];
L4: for($j=$im0+1;$j<$nmode;$j++) {
  for($ig=0;$ig<@ignored_modes;$ig++) {
    next L4 if($ignored_modes[$ig] == $j);
  }
  if($freq[$j] > $max_freq) {
    $max_freq = $freq[$j];
  } elsif($freq[$j] < $min_freq) {
    $min_freq = $freq[$j];
  }  
}
print "min_freq = $min_freq\n";
print "max_freq = $max_freq\n";

if($set_nrep eq 'no'){
  $nrep = @rep_names;
}
$multi_file = 'no';
$multi_file = 'yes' if(@rep_names > $nrep);

$level_length = 1.0/$nrep;

# gnuplot
$rest_rep = @rep_names;
$num = 0;
while($rest_rep > 0) {
$rest_rep -= $nrep;
$num ++;
if($multi_file eq 'no') {
  $epsf = "freq.eps"; 
} else {
  $epsf = "freq$num.eps"; 
}
$file_gpi = "$$.gpi";
open(OUT,">$file_gpi");
print OUT "set terminal postscript eps enhanced color solid 24\n";
print OUT "set output \"$epsf\"\n";
print OUT "set title 'Vibrational Analysis'\n";
print OUT "set size $width,$height\n";
print OUT "set noxtics\n";
print OUT "set ylabel \"Frequency (cm^{-1})\"\n";
print OUT "set nokey\n";
$label_pos = $max_freq*(1.0+0.1/$height); 
$istart = $nrep*($num-1);
$iend   = $nrep*$num;
$iend = @rep_names if($iend > @rep_names);
for($i=$istart;$i<$iend;$i++) {
  $ishift = $i - $nrep*($num-1);
  $start = $level_length*($ishift+0.2);
  $end   = $level_length*($ishift+0.7);
  $line_type = $i+1;
  for($j=0;$j<@{$freq_rep[$i]};$j++) {
    $no = $j+1;
    print OUT "set label '$no' at $start,$freq_rep[$i][$j] right\n";
    print OUT "set arrow  from $start,$freq_rep[$i][$j] to $end,$freq_rep[$i][$j] nohead lt $line_type  lw 4\n";
    $freq_int = int($freq_rep[$i][$j]+0.5);
    print OUT "set label '$freq_int' at $end,$freq_rep[$i][$j] left\n";
  }
  $act_names[$i] = 'IR\&R' if($act_names[$i] eq 'IR&R');
  print OUT "set label '$rep_names[$i] $act_names[$i]' at $start,$label_pos left\n";
}
$max_freq *= 1.0+0.2/$height;
print OUT "plot [0:1][0:$max_freq] 0";
close(OUT);

`gnuplot $file_gpi`;
`rm -f $file_gpi\n`;

}

exit;
