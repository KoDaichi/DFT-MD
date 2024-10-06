#!/usr/bin/perl
#
# A visualizaion program for molcular dynamics
# Version : 1.00
#
# Copyright (C) 2004, Takenori YAMAMOTO
# Institute of Industrial Science, The University of Tokyo
#

$bohr = 0.5291772480;
$eps  = 0.01;

if(@ARGV < 1 ) {
  print "Usage: dynm2tr2.pl nfdynm.data [ control.inp ]\n";
  print "Format of control.inp:\n";
  print " origin  X Y Z\n";
  print " vector1 A1x A1y A1z\n";
  print " vector2 A2x A2y A2z\n";
  print " vector3 A3x A3y A3z\n";
  exit;
}
$file = shift; 
$control_file_exists = 'no';
if(@ARGV > 0) {
  $control_file_exists = 'yes';
  $file_control = shift @ARGV;
}

@cell_origin = (0,0,0);
$vec1_exists = 'no';
$vec2_exists = 'no';
$vec3_exists = 'no';

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
  }
}
close(IN);
}

open(IN,$file);
<IN>;
for($i=0;$i<3;$i++) {
  ($dummy,$line) = split('=',<IN>);
  @{$vec[$i]} = split(' ',$line);
  if($vec1_exists eq 'no' && $vec2_exists eq 'no' && $vec3_exists eq 'no') {
    @{$cell_vector[$i]} = @{$vec[$i]};
  } 
  #print "vec = @{$vec[$i]}\n";
  #print "cell_vec = @{$cell_vector[$i]}\n";
}
@line = split(' ',<IN>);
$ntyp = $line[3];
$natm = $line[6];
print "ntyp= $ntyp, natm= $natm\n";
@type = ();
while(@type < $natm) {
  @line = split(' ',<IN>);
  shift @line;
  shift @line;
  push @type, @line;
}
for($i=0;$i<$ntyp;$i++) {
  @line = split(' ',<IN>);
  $name[$i] = $line[4];
}
<IN>;
open(OUT,">dynm.tr2");
$n = 1;
while($line = <IN>) {
  print "MD step $n\n";
  $num_atoms = 0;
  $text = '';
  for($i=0;$i<$natm;$i++) {
    @line = split(' ',<IN>);
    $name = $name[$type[$i]-1];
    for($j=0;$j<3;$j++) {
      $pos[$j] = $line[$j+1];
      $force[$j] = $line[$j+4];
    }
    #print "pos = @pos, force = @force\n";
    @duplpos = ();
    &duplication(\@pos,\@duplpos);
    $num_atoms += @duplpos;
    for($j=0;$j<@duplpos;$j++) {
      $pos[0] = $duplpos[$j][0]*$bohr;
      $pos[1] = $duplpos[$j][1]*$bohr;
      $pos[2] = $duplpos[$j][2]*$bohr;
      $text = $text . "$name @pos @force\n";
    }
  }
  print OUT "$num_atoms\n";
  print OUT "label = \"MD step $n\"\n";
  print OUT $text;
  $n++;
}
close(IN);
close(OUT);

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

sub duplication($$) {
  my ($pos,$duplpos) = @_;
  my ($n1,$n2,$n3,$k,@xyz,@red);

  my @cell_rvector =();
  &get_rvector(\@cell_vector,\@cell_rvector);
  my $n=0;
  for($n1=-10;$n1<=10;$n1++) {
  for($n2=-10;$n2<=10;$n2++) {
  for($n3=-10;$n3<=10;$n3++) {
    for($k=0;$k<3;$k++) {
      $xyz[$k] = $pos->[$k] + $vec[0][$k]*$n1 + $vec[1][$k]*$n2 + $vec[2][$k]*$n3;
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
      $duplpos->[$n][0] = $xyz[0];
      $duplpos->[$n][1] = $xyz[1];
      $duplpos->[$n][2] = $xyz[2];
      $n++;
    }
  }}}
}

sub get_rvector($$) {
  my ($vec,$rvec) = @_;
  my ($i,$i1,$i2,$i3,$j,$j1,$j2,$j3,$vol);
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
