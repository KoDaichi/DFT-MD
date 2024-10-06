#!/usr/bin/perl -w
#
# Density-of-states plot
# Version : 3.00
#
# Copyright (C) 2003-2004, Takenori YAMAMOTO
# Institute of Industrial Science, The University of Tokyo
# Copyright (C) 2005-2006, Takenori YAMAMOTO
# AdvaceSoft Corporation
#

use Math::Trig;

$file_dos = "dos_$$.dat";
$file_gpi = "dos_$$.gpi";

$set_erange = 'no';
$set_einc = 'no';
$set_dosrange = 'no';
$set_dosinc = 'no';
$set_fermi = 'no';
$set_color = 'no';
$set_epsf = 'yes';
$set_data = 'no';
$width=1;
$font_size=14;
$mode='total';
$title='';

if(@ARGV<1) {
  print "Version: 3.00\n";
  print "Usage: dos.pl DosData -erange=Emin,Emax -einc=dE -dosrange=DOSmin,DOSmax -dosinc=dDOS -title=STRING -with_fermi -width=SIZE -font=SIZE -color -mode={total|layer|atom|projected} -epsf={yes|no} -data={yes|no}\n";
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
    } elsif ($s =~/-dosrange/) {
      ($dummy,$s)=split('=',$s);
      ($dosmin,$dosmax)=split(',',$s);
      $set_dosrange = 'yes';
    } elsif ($s =~/-dosinc/) {
      ($dummy,$dosinc)=split('=',$s);
      $set_dosinc = 'yes';
    } elsif ($s =~/-title/) {
      ($dummy,$title)=split('=',$s);
    } elsif ($s =~/-with_fermi/) {
      $set_fermi = 'yes';
    } elsif($s =~/-width/) {
      ($dummy,$width)=split('=',$s);
    } elsif($s =~/-font/) {
      ($dummy,$font_size)=split('=',$s);
    } elsif ($s =~/-color/) {
      $set_color = 'yes';
    } elsif ($s =~/-mode/) {
      ($dummy,$mode)=split('=',$s);
    } elsif ($s =~/-epsf/) {
      ($dummy,$s)=split('=',$s);
      $set_epsf = &yes_or_no($s,$set_epsf);
    } elsif ($s =~/-data/) {
      ($dummy,$s)=split('=',$s);
      $set_data = &yes_or_no($s,$set_data);
    }
  }
}

sub yes_or_no($$) {
  my ($s,$defalt) = @_;
  if($s =~/^Y/ or $s =~/^y/) {
    return 'yes';
  } elsif($s =~/^N/ or $s =~/^n/) {
    return 'no';
  }
  return $defalt;
}

$file = $ARGV[0];
if($mode eq 'total') {
  &plot_dos($file,"density_of_states.eps",$title);
  &split_tdos_data($file) if($set_data eq 'yes');
} elsif($mode eq 'layer' or $mode eq 'atom' or $mode eq 'projected') {
  &plot_ldos($file,$mode);
} else {
  print "mode is invalid, which should be total, layer, or atom.";
}

exit;

# Subroutines
sub plot_dos($$$) {
my ($dos_file,$eps_file,$title) = @_;

return if($set_epsf eq 'no');

open(IN,$dos_file);
$line=<IN>;
$nspin = 1;
if($line=~/dos_up/) {
  $nspin = 2;
}
$ndos=0;
if($nspin==1) {
  while($line=<IN>) {
    last if($line=~/END/);
    @line = split(' ',$line);
    $energy[$ndos] = $line[3];
    $dos[$ndos] = $line[4];
    $sum[$ndos] = $line[5];
    $ndos++;
  }
} elsif($nspin == 2) {
  while($line=<IN>) {
    last if($line=~/END/);
    @line = split(' ',$line);
    $energy[$ndos] = $line[4];
    $dos_up[$ndos] = $line[5];
    $dos_down[$ndos] = $line[6];
    $sum_up[$ndos] = $line[7];
    $sum_down[$ndos] = $line[8];
    $sum_total[$ndos] = $line[9];
    $ndos++;
  }
}
close(IN);
 
# output DOS data
open(OUT,">$file_dos");
if($nspin == 1) {
  for($i=0;$i<$ndos;$i++) {
    print OUT "$energy[$i] $dos[$i] $sum[$i]\n";
  }
} elsif($nspin == 2) {
  for($i=0;$i<$ndos;$i++) {
    print OUT "$energy[$i] $dos_up[$i] $dos_down[$i] $sum_up[$i] $sum_down[$i] $sum_total[$i]\n";
  }
}
close(OUT);

# gnuplot
open(OUT,">$file_gpi");
if($set_color eq 'no') {
  print OUT "set terminal postscript eps enhanced solid $font_size\n";
} else {
  print OUT "set terminal postscript eps enhanced color solid $font_size\n";
}
print OUT "set output \"$eps_file\"\n";
print OUT "set size $width,1\n";
print OUT "set xlabel \"Energy (eV)\"\n";
print OUT "set ylabel \"DOS (states/eV)\"\n";
print OUT "set nokey\n";
if(length($title) > 0) {
  print OUT "set title \"$title\"\n";
}
if($set_fermi eq 'yes') {
  print OUT "set arrow  from 0, graph 0 to 0,graph 1 nohead lt 0\n";
}
if($set_einc eq 'yes') {
  print OUT "set xtics $einc\n";
}
if($set_dosinc eq 'yes') {
  print OUT "set ytics $dosinc\n";
}
print OUT "plot ";
if($set_erange eq 'yes') {
  print OUT "[$emin:$emax] ";
} else {
  print OUT "[:] ";
}
if($set_dosrange eq 'yes') {
  print OUT "[$dosmin:$dosmax] ";
} else {
  print OUT "[:] ";
}
print OUT "\\\n";
if($nspin == 1) {
  print OUT "\"$file_dos\" u 1:2 w l\n";
} elsif($nspin == 2) {
  print OUT "\"$file_dos\" u 1:2 w l lt 1, \"$file_dos\" u 1:(-(\$3)) w l lt 3, 0 lt -1\n";
}
close(OUT);

`gnuplot $file_gpi`;
`rm -f $file_dos\n`;
`rm -f $file_gpi\n`;
}

sub plot_ldos($$) {
  my ($file,$mode) = @_;
  my @files = ();
  my @info = ();
  &split_file($file,$mode,\@files,\@info);
  for(my $i=0;$i<@files;$i++) {
    my $eps_file = $files[$i];
    $eps_file =~s/data/eps/;
    my $title = $info[$i];
    my @line = split(' ',$title);
    if($mode eq 'layer') {
      $title = "Layer $line[3]";
    } elsif($mode eq 'atom') {
      $title = "Atom $line[3]";
    } elsif($mode eq 'projected') {
      $title = "PDOS atom $line[2] l=$line[4] m=$line[6] t=$line[8]";
    }
    &plot_dos($files[$i],$eps_file,$title);
    if($set_data eq 'no') {
      `rm $files[$i]`;
    }
  }
}

sub split_file($$$) {
  my ($file,$mode,$files,$info) = @_;
  my ($i,$tag,$line,$f,@line,$num,$l,$m,$t);

  if($mode eq 'layer') {
    $tag = 'LAYERDOS';
  } elsif($mode eq 'atom') {
    $tag = 'ALDOS';
  } elsif($mode eq 'projected') {
    $tag = 'PDOS';
  }

  open(IN,$file);
  $i=0;
  while($line = <IN>) {
    if($line =~/$tag/) {
      @line = split(' ',$line);
      if($mode eq 'layer') {
         $info->[$i] = $line;
	 $num = &id($line[3]);
         $f = "dos_l$num.data";
      } elsif($mode eq 'atom') {
         $info->[$i] = $line;
	 $num = &id($line[3]);
         $f = "dos_a$num.data";
      } elsif($mode eq 'projected') {
         $info->[$i] = $line;
	 $num = &id($line[2]);
	 $l = $line[4];
	 $m = &id2($line[6]);
	 $t = $line[8];
         $f = "dos_a${num}l${l}m${m}t${t}.data";
      }
      $files->[$i++] = $f;
      open(OUT,">$f");
      while($line=<IN>) {
        if($line =~/END/) {
          last;
	} else {
          print OUT $line;
	}
      }
      close(OUT);
    }
  }
  close(IN);
}

sub id($) {
  my ($i) = @_;
  my $digit = 3;
  my $len = length($i);
  my $num = $i;
  for(my $j=0;$j<$digit-$len;$j++) {
    $num = '0' . $num;
  }
  return $num;
}

sub id2($) {
  my ($i) = @_;
  my $digit = 2;
  my $len = length($i);
  my $num = $i;
  for(my $j=0;$j<$digit-$len;$j++) {
    $num = '0' . $num;
  }
  return $num;
}

sub split_tdos_data($) {
  my ($file) = @_;
  my ($line);

  open(IN,$file);
  open(OUT,">tdos.data");
  while($line=<IN>) {
    if($line =~/END/) {
      last;
    } else {
      print OUT $line;
    }
  }
  close(OUT);
}
