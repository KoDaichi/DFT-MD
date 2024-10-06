#!/usr/bin/perl -w

use strict;

our $width = 1.0;
our $keep_gpi;
our $fsize = 18;
our $help;
our @vlcr=();
our $vlcrfile;
our $wf;
our $efermi;
our $vmax = -10000000000;
our $mono;
our $vmaxpos;
our $nolabel;

our $nflat = 5;
our $cflat = 0.01;
our @flat;
our $flatregion_found=0;
our $center;
our $output = "workfunc.eps";

&parse_args;
&parse_wf_data;
&run_gnuplot;

sub print_usage {
    print "Usage : workfunc.pl nfvlcr_av.data [OPTIONS]\n";
}

sub parse_args {
    use Getopt::Std;
    use Getopt::Long;
    Getopt::Long::Configure('bundling');
    GetOptions(
    'o|output=s' => \$output,
    'width=f'=>\$width,
    'mono'=>\$mono,
    'nflat=i'=>\$nflat,
    'cflat=f'=>\$cflat,
    'gpi'=>\$keep_gpi,
    'nolabel'=>\$nolabel,
    'efermi=f' => \$efermi,
    );
    my $len=@ARGV;
    if($len==0 or defined($help)) {
        &print_usage;
        exit;
    }
    $vlcrfile=$ARGV[0];
    die("file : $vlcrfile does not exist \n") if(! -e $vlcrfile);
}

sub strip {
    my $str = $_[0];
    $str =~ s/^\s+//;
    $str =~ s/$str\s+$//;
    return $str;
}

sub parse_wf_data {
    open VLC,"<$vlcrfile";
    my @list = <VLC>;
    my $ndata = @list;
    $flatregion_found=0;
    for(my $i=0 ; $i<$ndata ; $i++) {
        my $line = $list[$i];
        $line = &strip($line);
        my @words=split(/\s+/,$line);
        if ($line =~ "#") {
            if ($line =~ "# Fermi" && !defined($efermi)){
                $efermi = $words[4];
            }
            next;
        }
        my @vl = ($words[0],$words[1]);
        push @vlcr, \@vl;
        $flat[$i] = 0;
    }

    for (my $i=0 ; $i<$ndata-$nflat ; $i++) {
       my $e0 = $vlcr[$i][1];
       my $test=1;
       for (my $j=1 ; $j<$nflat-1 ; $j++) {
           my $e1 = $vlcr[$i+$j][1];
           my $diff = abs(($e0-$e1)/$e0);
           if (!($diff<$cflat)){
               $test = 0;
               last;
           }
       }
       if ($test){
           $flat[$i] = 1;
           $flatregion_found = 1;
       }
    }

    if (!$flatregion_found) {
        print "failed to find a flat region in the local potential file : $vlcrfile\n";
    } else {
        my $nflat_region=0;
        $center = 0.0;
        my $nflat_max;
        for (my $i=0 ; $i<$ndata ; $i++) {
            next if(!$flat[$i]);
            if ($vmax<$vlcr[$i][1]) {
                $vmax = $vlcr[$i][1];
            }
            $nflat_region++;
            $center += $vlcr[$i][0];
            $nflat_max=$i;
        }
        $vmaxpos = $vlcr[$nflat_max][0];
        $wf = $vmax-$efermi;
        $center /= $nflat_region;
        print "estimated work function : $wf eV\n";
    }

}

sub run_gnuplot {
    my $gpi = "workfunc.gpi"; 
    open GP, ">$gpi";
    if(defined($mono)){
        print GP "set term post eps solid enhanced $fsize\n";
    } else {
        print GP "set term post eps color solid enhanced $fsize\n";
    }
    print GP "set output \"$output\"\n";
    print GP "set encoding iso\n";
    print GP "set xrange[0:]\n";
    print GP "set xlabel 'distance along the z-axis (\\305)'\n";
    print GP "set ylabel 'energy (eV)\n";
    if ($flatregion_found && !defined($nolabel)) {
        my $epos = $efermi+abs($efermi-$vmax)*0.5;
        my $wfstr = sprintf ("%.2f",$wf);
        my $labelpos = $center+$vmaxpos*0.01;
        print GP "set arrow from first $center,first $efermi to first $center,first $vmax heads lt -1 \n";
        print GP "set label '{/Symbol f} = $wfstr eV' at $labelpos,$epos\n";
    }
    print GP "set size $width,1\n";
    print GP "plot '$vlcrfile' using 1:2 w l not , $efermi w l not lt 0\n";
    close(GP);

    my $osname = $^O;
    my $gnuplot="";
    my $rm = "";
    if ($osname =~ "lin"){ 
        $gnuplot="gnuplot";
        $rm = "rm -f";
    } elsif ($osname =~ "MS"){
        $gnuplot="wgnuplot.exe";
        $rm = "del";
    }
    `$gnuplot $gpi\n`;
    if (!defined($keep_gpi)) {
        `$rm $gpi\n`;
    }
}

