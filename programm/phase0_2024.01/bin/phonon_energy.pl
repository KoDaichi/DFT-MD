#!/usr/bin/perl -w

use strict;

our $ha2ev = 27.2113834;
our $kb_ev = 8.617343e-5;
our $kb_ha = $kb_ev/$ha2ev;

our $modefile;
our $tmpenergyfile = "phonon_energy.data";
our $Tmax = 3000;
our $Tmin = 0;
our $nT = 100;

our $width=1.0;
our $trange="[0:3000]";
our $tinc;
our $cinc;
our $einc;

our $help;
our $mono;
our $fsize=18;
our $keep_gpi;

our %band_data;
our $nstates;
our $nqvec;
our $natm;

&parse_args;
&parse_band_data;
&write_data_file;
&run_gunplot;

exit;

sub parse_args {
    use Getopt::Std;
    use Getopt::Long;
    Getopt::Long::Configure('bundling');
    GetOptions(
    'w|width=f' => \$width,
    't|trange=s' => \$trange,
    'tinc=f'=> \$tinc,
    'n|nT=i'=> \$nT,
    'einc=f'=> \$einc,
    'cinc=f'=> \$cinc,
    'f|font=i'=> \$fsize,
    'mono'=> \$mono,
    'h|help'=> \$help,
    'gpi'=> \$keep_gpi,
    );
    my $len=@ARGV;
    if ($len==0 or defined($help)) {
        &print_usage;
        exit;
    }
    $modefile = $ARGV[0];
    if (! -e $modefile ){
        die("file : $modefile does not exist\n");
    }
    if (defined($trange)) {
        if ($trange !~ m/(\[-*\+*\d*:-*\+*\d*\])/){
            die("invalid --trange option : $trange\n");
        }
        if ($trange eq "[:]") {
            $trange = "[0:3000]"
        } else {
            my $tmpstr = substr($trange,1,length($trange)-2);
            my @tr = split(/:/,$tmpstr);
            $Tmin = $tr[0];
            $Tmax = $tr[1];
        }
    }
}

sub print_usage {
    print "Usage : phonon_energy.pl mode.data [OPTIONS]\n";
    print "\n";
    print "OPTIONS\n";
    print "--width  or -w : specify the width of the graph. defaults to 1. \n";
    print "--trange or -t : specify the range of the temperatur by the form [tmin:tmax].\n";
    print "--nT     or -n : specify the number of temperature points.\n";
    print "--font   or -f : specify the fontsize to be used in the graph. defaults to 18.\n";
    print "--mono   or -m : specify this option in order to create monochrome figures.\n";
    print "--tinc         : specify the tic for the x-axis.\n";
    print "--einc         : specify the tic for the y-axis (energy).\n";
    print "--cinc         : specify the tic for the y-axis (specific heat).\n";
}

sub strip {
    my $str = $_[0];
    $str =~ s/^\s+//;
    $str =~ s/$str\s+$//;
    return $str;
}

sub parse_band_data {
    open MD,"<$modefile" or die("$!");
    my @list=<MD>;
    my @eigs = ();
    my $qid=0;
    my $currmode=0;
    my $header_resolved=0;
    my $weight=1;
    foreach my $line (@list) {
        $line = &strip($line);
        if ($line =~ "---") {
            next;
        }
        my @words=split(/\s+/,$line);
        if ($line =~ "Nmode") {
            if(scalar(@words)<6) {
                die("invalid line : $line read from file : $modefile\n");
            }
            $nstates = $words[1];
            $natm = $words[3];
            $nqvec = $words[5];
            $header_resolved = 1;
            next;
        }
        if (!$header_resolved) {
            next;
        }
        if ($line =~ "iq=") {
           if(scalar(@words)<10) {
                die("invalid line : $line read from file : $modefile\n");
           }
           $qid = $words[1];
           die("weight undefined for q-point no. $qid") if (scalar(@words)<11);
           $weight = $words[10];
           next;
       }
       if ($line =~ "hbarW") {
           if (scalar(@words)<9){
                die("invalid line : $line read from file : $modefile\n");
           }
           my $e = $words[1];
           push(@eigs,$e);
           $currmode++;
           if($currmode == $nstates) {
               my @band=($weight);
               foreach my $eig (@eigs) {
                   push(@band,$eig);
               }
               $band_data{$qid} = \@band;
               $currmode=0;
               @eigs=();
               next;
           }
       }
    }
    close(MD);
}

sub run_gunplot {
    my $gpi = "phonon_energy.gpi";
    open GP,">$gpi" or die("$!");

    if (defined($mono)) {
        print GP "set term post eps enhanced $fsize\n";
    } else {
        print GP "set term post eps color solid enhanced $fsize\n";
    }

    print GP "set output \"phonon_energy.eps\"\n";

    if (defined($tinc)) {
        print GP "set xtic $tinc\n";
    }

    if (defined($einc)) {
        print GP "set ytic $einc\n";
    }

    print GP "set xrange $trange \n";

    print GP "set size $width,1\n";
    print GP "set xlabel \"temperature (K)\"\n";
    print GP "set ylabel \"energy (eV)\"\n";
    print GP "set y2label \"entropy (eV/K)\"\n";
    print GP "set ytics nomirror\n";
    print GP "set y2tic\n";
    print GP "plot \"$tmpenergyfile\" using 1:2 w l t \"Internal Energy\",";
    print GP "\"$tmpenergyfile\" using 1:3 w l t \"Free Energy\",";
    print GP "\"$tmpenergyfile\" using 1:4 axis x1y2 w l t \"Entropy\"\n";
    print GP "\n";
    close GP;

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

    $gpi = "phonon_Cv.gpi";
    open GP,">$gpi" or die("$!");

    if (defined($mono)) {
        print GP "set term post eps enhanced $fsize\n";
    } else {
        print GP "set term post eps color solid enhanced $fsize\n";
    }

    print GP "set output \"phonon_Cv.eps\"\n";

    if (defined($tinc)) {
        print GP "set xtic $tinc\n";
    }

    if (defined($cinc)) {
        print GP "set ytic $cinc\n";
    }

    print GP "set xrange $trange \n";

    print GP "set size $width,1\n";
    print GP "set xlabel \"temperature (K)\"\n";
    print GP "set ylabel \"C_{v} (k_{B}/atom)\"\n";
    print GP "plot \"$tmpenergyfile\" using 1:5 w l not ";
    close GP;

    `$gnuplot $gpi\n`;

    if (!defined($keep_gpi)) {
        `$rm $gpi\n`;
    }

}

sub write_data_file {
    open EN, ">$tmpenergyfile" or die("$!");
    print EN "# T (K) Internal Energy (eV) Free energy (eV) Entropy (eV/K) Cv (kB/atom)\n";
    my $dT = ($Tmax - $Tmin)/$nT;
    for (my $t=0 ; $t<=$nT ; $t++) {
        my $temperature = $Tmin + $dT * $t;
        my @edata = &get_phonon_energy_etc($temperature);
        print EN "$temperature $edata[0] $edata[1] $edata[2] $edata[3] \n";
    }
    close(EN);
}

sub get_phonon_energy_etc {
    my $temperature = $_[0];
    my $u = 0.0;
    my $f = 0.0;
    my $c = 0.0;
    for (my $q=1 ; $q<=$nqvec ; $q++) {
        my $band = $band_data{$q};
        my $weight = $band->[0];
        for (my $ie=0; $ie<$nstates ; $ie++) {
            my $e = $band->[$ie+1];
            next if($e<=0);
            if($temperature>0.1) {
                my $x = $e/($kb_ha*$temperature);
                $u += $weight*(0.5*$e+$e/(exp($x)-1.0));
                $f += $weight*(0.5*$e+$kb_ha*$temperature*log(1-exp(-$x)));
                $c += $weight*$x**2*exp($x)/(exp($x)-1)**2;
            } else {
                $u += $weight * $e * 0.5;
                $f += $weight * $e * 0.5;
                $c += 0.0;
            }
        }
    }
    $u *= $ha2ev;
    $f *= $ha2ev;
    my $s=0.0;
    $s = ($u-$f)/$temperature if ($temperature>0.1);
    $c /= $natm;
    return ($u,$f,$s,$c);
}

