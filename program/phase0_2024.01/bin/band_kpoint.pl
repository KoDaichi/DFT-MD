#!/usr/bin/perl -w
#
# k-point generator for energy band plot
# version 2.0
#
# input file:
# distance
# b1x b2x b3x
# b1y b2y b3y
# b1z b2z b3z
# k1 k2 k3 denom # X
# k1 k2 k3 denom # Gamma
# k1 k2 k3 denom # L
# ....
#
#  Contact address :  Phase System Consortium
#

my @dummy;

$dummy2  = "";
$outfile = 'none';
$zdiv = 0;
$zshift = 0;

if (@ARGV<1) {
    print "Usage: \n";
    print "Usage: band_kpoint.pl InputFile -outfile=OUTFILE -zdiv=N1 -zshift=N2 \n";
    print "Format of input file: \n";
    print "distance \n";
    print "b1x b2x b3x \n";
    print "b1y b2y b3y \n";
    print "b1z b2z b3z \n";
    print "k1 k2 k3 denom # X \n";
    print "k1 k2 k3 denom # Gamma \n";
    print "k1 k2 k3 denom # L \n";
    print ".... \n";
    exit;
} else {
    foreach $s ( @ARGV ) {
	if ($s =~/-outfile/) {
	    ($dummy2,$outfile)=split('=',$s);
        } elsif ($s =~/-zdiv/) {
	    ($dummy2,$zdiv)=split('=',$s);
        } elsif ($s =~/-zshift/) {
	    ($dummy2,$zshift)=split('=',$s);
        }
    }
}

if ( $outfile eq 'none' ) {
    $outfile = "kpoint.data"
}

open(IN,$ARGV[0]);
$distance = <IN>;
chomp($distance);
@{$rlvec[0]} = split(' ',<IN>);
@{$rlvec[1]} = split(' ',<IN>);
@{$rlvec[2]} = split(' ',<IN>);

$i=0;
while($line = <IN>) {
  ($line,@dummy) = split('#',$line);
  @{$skp[$i]} = split(' ',$line);
  die "bad input: i=$i skp=@{$skp[$i]}\n" if(@{$skp[$i]} != 4);
  for($j=0;$j<3;$j++) {
    $skv[$i][$j] = $skp[$i][$j]/$skp[$i][3];
  }
  $i++;
}
close(IN);
$nsk = @skv;

# Division numbers will be calculated by the inputed distance.
for($i=1;$i<$nsk;$i++) {
  @dkc=(0,0,0);
  for($j=0;$j<3;$j++) {
    for($k=0;$k<3;$k++) {
      $dkc[$j] += $rlvec[$j][$k]*($skv[$i][$k]-$skv[$i-1][$k]);
    }
  }
  $dst = sqrt($dkc[0]**2+$dkc[1]**2+$dkc[2]**2);
  print "Distance of $i = $dst\n";
  $ndiv[$i-1] = int($dst/$distance);
}

print "division numbers = @ndiv\n";

# orientation vectors
for($i=0;$i<$nsk-1;$i++) {
  for($j=0;$j<3;$j++) {
    $orivec[$i][$j] = ($skv[$i+1][$j]-$skv[$i][$j])/$ndiv[$i];
  }
  print "i=$i @{$orivec[$i]}\n";
}

# generate k-points
$nkv=0;
$deno[$nkv] = $skp[0][3]*$skp[1][3]*$ndiv[0];
for($k=0;$k<3;$k++) {
  $kpt[$k] = $skv[0][$k];
  if($kpt[$k] > 0) {
    $nume[$nkv][$k] = int($kpt[$k]*$deno[$nkv]+0.5);
  } else {
    $nume[$nkv][$k] = int($kpt[$k]*$deno[$nkv]-0.5);
  }
}
$n=$nkv+1;
print "$n : @kpt\n    => $nume[$nkv][0]/$deno[$nkv] $nume[$nkv][1]/$deno[$nkv] $nume[$nkv][2]/$deno[$nkv]\n";
$nkv++;
for($i=1;$i<$nsk;$i++) {
  @kpt = @{$skv[$i-1]};
  for($j=0;$j<$ndiv[$i-1];$j++) {
    $deno[$nkv] = $skp[$i-1][3]*$skp[$i][3]*$ndiv[$i-1];
    for($k=0;$k<3;$k++) {
      $kpt[$k] += $orivec[$i-1][$k];
      if($kpt[$k] > 0) {
        $nume[$nkv][$k] = int($kpt[$k]*$deno[$nkv]+0.5);
      } else {
        $nume[$nkv][$k] = int($kpt[$k]*$deno[$nkv]-0.5);
      }
    }
    $n=$nkv+1;
    print "$n : @kpt\n    => $nume[$nkv][0]/$deno[$nkv] $nume[$nkv][1]/$deno[$nkv] $nume[$nkv][2]/$deno[$nkv]\n";
    $nkv++;
  }
}

# write data in PHASE format
open(OUT,">$outfile");

if ( $zdiv <= 1 ) {
    print OUT "$nkv $nkv\n";
    for($i=0;$i<$nkv;$i++) {
	print OUT "$nume[$i][0] $nume[$i][1] $nume[$i][2] $deno[$i] 1\n";
    }
} else {
    $nkv2 = $nkv * $zdiv;
    print OUT "$nkv2 $nkv2 \n";

    if ( $zdiv%2 == 0 ) {
	$dz = 1.0 /2.0;
    } else {
	$dz = ($zdiv -1) /$zdiv /2.0;
    }

    $kzshift = 0.0;
    if ( $zshift > 0 ) {
	$kzshift = 0.5;
    }
    for($i=0;$i<$nkv;$i++) {
	$kx = $nume[$i][0] /$deno[$i];
	$ky = $nume[$i][1] /$deno[$i];
	$kz = $nume[$i][2] /$deno[$i];
	for($j=0;$j<$zdiv;$j++) {
	    $kz2 = $kz +( $j +$kzshift )/$zdiv -$dz;
	    printf OUT ("%25.20f%25.20f%25.20f\n", $kx, $ky, $kz2);
	}
    }

}
close(OUT);
