#!/usr/bin/perl -w

use strict;

our $ha2ev      = 27.21139615;
our $ha2cminv  = $ha2ev/1.23984185e-4;
our $ha2thz    = $ha2ev/1.23984185e-4/33;
our %unitfactor = (
"cm-1"=>$ha2cminv,
"mHa"=>1000.0,
"meV"=>$ha2ev*1000,
"THz"=>$ha2thz
);

our $modefile;
our $tmpbandfile = "phonon_band.data";
our $controlfile="bandkpt.in";
our $title;
our $width=1.0;
our $units="cm-1";
our $erange;
our $einc;
our $fsize=18;
our $keep_data;
our $keep_gpi;
our $mono;
our $help;
our $ptype = "line";

our $emin=1e30;
our $emax=-1e30;
our $nstates=0;
our $qidmax=0;

our %band_data;

our $qmax;
our %qlen_at_sp;

our @special_kpoints;
our $nspk=-1;

&parse_args;
&parse_control;
&parse_band_data;
&write_data_file;
&run_gnuplot;

exit;

sub parse_args {
    use Getopt::Std;
    use Getopt::Long;
    Getopt::Long::Configure('bundling');
    GetOptions(
    't|title=s' => \$title,
    'w|width=f' => \$width,
    'c|control=s' => \$controlfile,
    'u|units=s' => \$units,
    'e|erange=s' => \$erange,
    'einc=f'=> \$einc,
    'k|keep'=> \$keep_data,
    'f|font=i'=> \$fsize,
    'p|ptype=s'=> \$ptype,
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

    if (defined($erange)) {
        if ($erange !~ m/(\[-*\+*\d*:-*\+*\d*\])/){
            die("invalid --erange option : $erange\n");
        }
    }

    if (!exists($unitfactor{$units})){
        die("option --units must be one of : mHa, meV, THz, cm-1\n");
    }

    if (!($ptype eq "line" or $ptype eq "circle")) {
        die("option --ptype must be one of : line or circle");
    }
}

sub parse_control {
    if(! -e $controlfile ) {
        print "WARN : file '$controlfile' does not exist.\n";
        print "WARN : special k-points will not be taken into account.\n";
        return;
    }
    open CTRL, "<$controlfile" or die("$!");
    my $distance = <CTRL>;
    chomp($distance);
    my @rlvec;
    @{$rlvec[0]} = split(' ',<CTRL>);
    @{$rlvec[1]} = split(' ',<CTRL>);
    @{$rlvec[2]} = split(' ',<CTRL>);
    while(my $line=<CTRL>) {
        my @words = split('#',$line);
        my $symbol = '';
        if (scalar(@words)>1) {
            chomp($words[1]);
            $symbol = $words[1];
        }
        my @kp = split(' ',$words[0]);
        if (scalar(@kp)>=4) {
            my $n1 = $kp[0];
            my $n2 = $kp[1];
            my $n3 = $kp[2];
            my $nd = $kp[3];
            $n1 /= $nd;
            $n2 /= $nd;
            $n3 /= $nd;
            my @special_k = ($symbol,$n1,$n2,$n3);
            push @special_kpoints,\@special_k;
        }
    }
    $nspk = @special_kpoints;

    close(CTRL);
}

sub get_symbol {
    my $symb = $_[0];
    if ($symb =~ "G") {
        return "{/Symbol G}";
    }
    return $symb;
}

sub parse_band_data {
    open MD,"<$modefile" or die("$!");
    my @list=<MD>;
    my @eigs = ();
    my $qid=0;
    my $nmode=0;
    my $currmode=0;
    my $header_resolved=0;
    my @qnow=(0,0,0);
    my @qprev=(0,0,0);
    my $dq=0.0;
    my $qlen=0.0;
    my $firstq=1;
    my @qfract = (0,0,0);
    my $eps = 1.0e-8;
    my $isp = 0;
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
            $nmode = $words[1];
            $nstates = $nmode;
            $header_resolved = 1;
            next;
        }
        if (!$header_resolved) {
            next;
        }
        if ($line =~ "iq=") {
           die("invalid line : $line read from file : $modefile\n") if (scalar(@words)<10) ;
           $qid = $words[1];
           chop($words[3]);
           chop($words[4]);
           chop($words[5]);
           chop($words[7]);
           chop($words[8]);
           chop($words[9]);
           $qfract[0] = $words[3];
           $qfract[1] = $words[4];
           $qfract[2] = $words[5];
           $qnow[0]=$words[7];
           $qnow[1]=$words[8];
           $qnow[2]=$words[9];
           if (!$firstq) {
               $dq = sqrt(($qnow[0]-$qprev[0])**2+($qnow[1]-$qprev[1])**2+($qnow[2]-$qprev[2])**2);
           } else {
               $dq = 0.0;
               $firstq=0;
           }
           $qlen += $dq;
           $qprev[0] = $qnow[0];
           $qprev[1] = $qnow[1];
           $qprev[2] = $qnow[2];
           for (my $i=0 ; $i<$nspk ; $i++) {
               if(abs($qfract[0]-$special_kpoints[$i][1])<$eps &&
                  abs($qfract[1]-$special_kpoints[$i][2])<$eps &&
                  abs($qfract[2]-$special_kpoints[$i][3])<$eps) {
                   $qlen_at_sp{$isp} = $qlen;
                   $isp++;
                   last;
               }
           }
           next;
       }
       if ($line =~ "hbarW") {
           if (scalar(@words)<9){
                die("invalid line : $line read from file : $modefile\n");
           }
           my $e = $words[1];
           push(@eigs,$e*$unitfactor{$units});
           $currmode++;
           if($currmode == $nmode) {
               my @band=($qlen);
               foreach my $eig (@eigs) {
                   push(@band,$eig);
                   if ($eig<$emin) {
                      $emin = $eig;
                   } elsif ($eig>$emax) {
                      $emax = $eig;
                   }
               }
               $band_data{$qid} = \@band;
               $currmode=0;
               @eigs=();
               next;
           }
       }
    }
    close(MD);
    $qidmax = $qid;
    $qmax = $qlen;
}

sub is_special_kpoint {
    my $eps = 1.0e-4;
    my @qv0 = $_[0];
    my @qv1 = $_[1];
    return abs(($qv0[0]-$qv1[1])**2+($qv0[1]-$qv1[1])**2+($qv0[2]-$qv1[2])**2)<$eps;
}

sub write_data_file {
    open BAND,">$tmpbandfile" or die("$!");
    for (my $iq=1 ; $iq<=$qidmax ; $iq++) {
        my $bdata = $band_data{$iq};
        my $qlen = $bdata->[0];
        print BAND "$qlen ";
        for (my $j=1 ; $j<$nstates+1 ; $j++){
            my $e = $bdata->[$j];
            print BAND "$e ";
        }
        print BAND "\n";
    }
    close(BAND);
}

sub run_gnuplot {
    my $gpi = "phonon_band.gpi";
    open GP,">$gpi" or die("$!");

    if (defined($mono)) {
        print GP "set term post eps solid enhanced $fsize\n";
    } else {
        print GP "set term post eps color solid enhanced $fsize\n";
    }
    print GP "set output \"phonon_band.eps\"\n";
    print GP "unset xtic\n";
    print GP "set size $width,1\n";

    for (my $isp=0 ; $isp<$nspk ; $isp++) {
        print GP "set arrow from $qlen_at_sp{$isp},graph 0 to $qlen_at_sp{$isp},graph 1 nohead lt -1\n";
    }

    if($nspk>0){
        print GP "set xtics (";
        for (my $isp=0 ; $isp<$nspk ; $isp++) {
            my $sp = $special_kpoints[$isp];
            my $lab1 = $sp->[0];
            print GP "\"$lab1\" $qlen_at_sp{$isp}";
            print GP ", " if ($isp!=$nspk-1);
        }
        print GP ")\n";
    }

    if (defined($einc)) {
        print GP "set ytic $einc\n";
    }

    if($units eq 'meV') {
        print GP "set ylabel \"Frequency (meV)\"\n";
    } elsif($units eq 'cm-1') {
        print GP "set ylabel \"Frequency (cm^{-1})\"\n";
    } elsif($units eq 'THz') {
        print GP "set ylabel \"Frequency (THz)\"\n";
    } else {
        print GP "set ylabel \"Frequency (mHa)\"\n";
    }

    if (defined($title)) {
        print GP "set title \"$title\"\n";
    }

    if(defined($erange)){
        print GP "plot [0\:$qmax]$erange ";
    } else {
        print GP "plot [0\:$qmax][$emin\:$emax*1.02] ";
    }

    my $i=0;
    for ($i=0 ; $i<$nstates ; $i++) {
        my $ndat=$i+2;
        my $pt = "w l lt 1";
        if ($ptype eq "circle"){
            $pt = "w p pt 7 lt 1";
        }
        if ($i==$nstates-1) {
            print GP "\"$tmpbandfile\" u 1:$ndat $pt not ";
        } else {
            print GP "\"$tmpbandfile\" u 1:$ndat $pt not, ";
        }
    }
    print GP "\n";
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

    if (!defined($keep_data)) {
        `$rm $tmpbandfile\n`;
    }
}

sub strip {
    my $str = $_[0];
    $str =~ s/^\s+//;
    $str =~ s/$str\s+$//;
    return $str;
}

sub print_usage {
    print "Usage : phonon_band.pl mode.data [OPTIONS]\n";
    print "\n";
    print "OPTIONS\n";
    print "--units   or -u : specify the unit of energy. One of : mHa, meV, THz or cm-1.\n";
    print "--ptype   or -p : specify the 'plot type'. either line or circle, defaults to line.\n";
    print "--width   or -w : specify the width of the graph. defaults to 1. \n";
    print "--control or -c : specify the control file for the special k-points. defaults to 'bandkpt.in'. \n";
    print "--erange  or -e : specify the range of the energy by the form [emin:emax].\n";
    print "--title   or -t : specify the title of the graph. \n";
    print "--font    or -f : specify the fontsize to be used in the graph. defaults to 18.\n";
    print "--keep    or -k : specify this option in order to keep the intermediate data file.\n";
    print "--mono    or -m : specify this option in order to create monochrome figures.\n";
    print "--einc          : specify the tic for the y-axis.\n";
}

