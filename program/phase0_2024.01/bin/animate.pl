#!/usr/bin/perl -w
#
# A visualizaion program for lattice vibration
# Version : 1.00
#
# Copyright (C) 2004, Takenori YAMAMOTO
# Institute of Industrial Science, The University of Tokyo
#

use Math::Trig;

if(@ARGV < 1 || @ARGV > 2) {
  print "Usage: animate.pl modes.data [ control.inp ]\n";
  print "Format of control.inp:\n";
  print " origin  X Y Z\n";
  print " vector1 A1x A1y A1z\n";
  print " vector2 A2x A2y A2z\n";
  print " vector3 A3x A3y A3z\n";
  print " max_displacement Umax\n";
  exit;
}

$file_modes = shift @ARGV;
$control_file_exists = 'no';
if(@ARGV > 0) {
  $control_file_exists = 'yes';
  $file_control = shift @ARGV;
}

$bohr = 0.5291772480;

$umax = 0.5;
$ndiv = 20;
$cycle = 1;
$eps = 1.e-3;

@cell_origin = (0,0,0);
$vec1_exists = 'no';
$vec2_exists = 'no';
$vec3_exists = 'no';

$dielectric_mode = 'no';

# read modes.data
open(IN,"$file_modes");
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
  <IN>;
  ($dummy,$omega[$j],$dummy,$dummy,$omega_ev[$j],$dummy,$dummy,$freq[$j]) = split(' ',<IN>);
  $omega[$j] =~ s/D/E/;
  $omega_ev[$j] =~ s/D/E/;
  $freq[$j] =~ s/D/E/;
  for($i=0;$i<$natom;$i++) {
    ($no,@{$u[$j][$i]}) = split(' ',<IN>);
    $norm[$j][$i] = 0;
    for($k=0;$k<3;$k++) {
      $u[$j][$i][$k] /= sqrt($mass[$i]);
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

# read control.inp
if($control_file_exists eq 'yes') {
open(IN,"$file_control");
while($line = <IN>) {
  ($type,@data) = split(' ',$line);
  if($type =~/origin/) {
    @cell_origin = @data; 
    #$origin_exists = 'yes';
  } elsif($type =~/vector1/) {
    @{$cell_vector[0]} = @data; 
    $vec1_exists = 'yes';
  } elsif($type =~/vector2/) {
    @{$cell_vector[1]} = @data; 
    $vec2_exists = 'yes';
  } elsif($type =~/vector3/) {
    @{$cell_vector[2]} = @data; 
    $vec3_exists = 'yes';
  } elsif($type =~/max_displacement/) {
    ($umax) = @data;
  }
}
close(IN);
}

if($vec1_exists eq 'no' && $vec2_exists eq 'no' && $vec3_exists eq 'no') {
  @cell_vector = @a;
}

# duplication of atoms in the primitive cell
@cell_rvector =();
&get_rvector(\@cell_vector,\@cell_rvector);
for($i=0;$i<$natom;$i++) {
  $count=0;
  for($n1=-10;$n1<=10;$n1++) {
  for($n2=-10;$n2<=10;$n2++) {
  for($n3=-10;$n3<=10;$n3++) {
    for($k=0;$k<3;$k++) {
      $xyz[$k] = $pos[$i][$k] + $a[0][$k]*$n1 + $a[1][$k]*$n2 + $a[2][$k]*$n3;
      $xyz[$k] -= $cell_origin[$k];
    }
    for($k=0;$k<3;$k++) {
      $red[$k] = $cell_rvector[$k][0]*$xyz[0] + 
                 $cell_rvector[$k][1]*$xyz[1] + 
                 $cell_rvector[$k][2]*$xyz[2];
    }
    if(-$eps < $red[0] && $red[0] < 1.0 + $eps &&
       -$eps < $red[1] && $red[1] < 1.0 + $eps &&
       -$eps < $red[2] && $red[2] < 1.0 + $eps ) {
      $dupl[$i][$count++] = [$n1,$n2,$n3];
    }
  }}}
}

# normalization of displacement vectors
print "umax = $umax\n";
print "max_norm = $max_norm\n";
$ratio = $umax/$max_norm;
for($j=0;$j<$nmode;$j++) {
  for($i=0;$i<$natom;$i++) {
    for($k=0;$k<3;$k++) {
      $u[$j][$i][$k] *= $ratio;
    }
  }
}

# write TRJ files    
$total = $ndiv*$cycle;
$natom_xyz = 0;
for($i=0;$i<$natom;$i++) {
  $natom_xyz += @{$dupl[$i]};
}
for($j=0;$j<$nmode;$j++) {
  $no = $j+1;
  open(OUT,">mode_$no" . ".tr2");
  for($t=0;$t<$total;$t++) {
    $step = $t+1;
    print OUT "$natom_xyz\n";
    printf OUT "label = \"step %5i/%5i\" omega(Ha)=%10.3e omega(eV)=%10.3e nu(cm-1)=%10.3f\n", $step, $total, $omega[$j], $omega_ev[$j], $freq[$j];
    $d = sin($t/$ndiv*2.0*pi);
    for($i=0;$i<$natom;$i++) {
      for($k=0;$k<3;$k++) {
        $xyz[$k] = $pos[$i][$k] + $u[$j][$i][$k]*$d;
      }
      for($n=0;$n<@{$dupl[$i]};$n++) {
        for($k=0;$k<3;$k++) {
          $xyz2[$k] = $xyz[$k]
                    + $a[0][$k]*$dupl[$i][$n][0]
                    + $a[1][$k]*$dupl[$i][$n][1]
                    + $a[2][$k]*$dupl[$i][$n][2] - $cell_origin[$k];
          $xyz2[$k] *= $bohr;
          $vec2[$k] = $u[$j][$i][$k]*$bohr;
        }
        printf OUT "%4s %22.14f %22.14f %22.14f %22.14f %22.14f %22.14f\n", $name[$i], @xyz2, @vec2;
      }
    }
  }
}

# write grid mol2 file
for($i=0;$i<3;$i++) {
for($k=0;$k<3;$k++) {
  $cell_vector[$i][$k] *= $bohr;
}}
open(OUT,">grid.mol2");
print OUT <<here;
@<TRIPOS>MOLECULE
crystalline unit cell
8 12 0
0 0 0 0


grid file
@<TRIPOS>ATOM
here
$nump = 0;
for($i1=0;$i1<=1;$i1++) {
  for($i2=0;$i2<=1;$i2++) {
    for($i3=0;$i3<=1;$i3++) {
      $nump++;
      for($k=0;$k<3;$k++) {
        $grid_point[$k] = $i1*$cell_vector[0][$k] + $i2*$cell_vector[1][$k] + $i3*$cell_vector[2][$k];
      }
      print OUT "$nump N @grid_point N.4 1 GLY 0.0000\n"
    }
  }
}
print OUT <<here;
@<TRIPOS>BOND
1 5 1
2 6 5
3 5 7
4 3 1
5 3 4
6 3 7
7 2 1
8 2 4
9 2 6
10 8 4
11 8 6
12 8 7
here
close(OUT);

#
sub get_rvector($$) {
  my ($vec,$rvec) = @_;
  for($i=0;$i<3;$i++) {
    $i1 = ($i+1)%3;
    $i2 = ($i+2)%3;
    for($j=0;$j<3;$j++) {
      $j1 = ($j+1)%3;
      $j2 = ($j+2)%3;
      $rvec->[$i][$j] = $vec->[$i1][$j1]*$vec->[$i2][$j2]-$vec->[$i1][$j2]*$vec->[$i2][$j1];
    }
  }
  $vol = $rvec->[0][0]*$vec->[0][0] + $rvec->[0][1]*$vec->[0][1] + $rvec->[0][2]*$vec->[0][2];
  for($i=0;$i<3;$i++) {
    for($j=0;$j<3;$j++) {
      $rvec->[$i][$j] /= $vol;
    }
  }
}
