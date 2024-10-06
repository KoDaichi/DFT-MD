#!/usr/bin/perl -w

use strict;

our $ha2ev      = 27.2113834;
our $mha2cminv  = $ha2ev/1000/1.23984185e-4;
our $mha2thz    = $ha2ev/1000/1.23984185e-4/33;
our %unitfactor = (
"cm-1"=>$mha2cminv,
"mHa"=>1.0,
"meV"=>$ha2ev,
"THz"=>$mha2thz
);

our $dosfile;
our $title;
our $width=1.0;
our $units="cm-1";
our $erange;
our $drange;
our $fsize = 18;
our $dinc;
our $einc;
our $keep_data;
our $keep_gpi;
our $mono;
our $help;

our $tmpdosfile = "phonon_dos.data";
our @dosdata=();

&parse_args;
&parse_dos_data;
&write_dos_data;
&run_gnuplot;

exit;

sub parse_args {
    use Getopt::Std;
    use Getopt::Long;
    Getopt::Long::Configure('bundling');
    GetOptions(
    't|title=s' => \$title,
    'w|width=f' => \$width,
    'u|units=s' => \$units,
    'e|erange=s' => \$erange,
    'd|drange=s' => \$drange,
    'dinc=f'=> \$dinc,
    'einc=f'=> \$einc,
    'k|keep'=> \$keep_data,
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

    $dosfile=$ARGV[0];
    if (! -e $dosfile ){
      die("file : $dosfile does not exist\n");
    }

    if (defined($erange)) {
        if ($erange !~ m/(\[-*\+*\d*:-*\+*\d*\])/){
            die("invalid --erange option : $erange\n");
        }
    }

    if (defined($drange)) {
        if ($drange !~ m/(\[-*\+*\d*:-*\+*\d*\])/){
            die("invalid --drange option : $drange\n");
        }
    }

    if (!exists($unitfactor{$units})){
        die("option --units must be one of : mHa, meV, THz, cm-1\n");
    }
}

sub strip {
    my $str = $_[0];
    $str =~ s/^\s+//;
    $str =~ s/$str\s+$//;
    return $str;
}

sub print_usage {
    print "Usage : phonon_dos.pl phdos.data [OPTIONS]\n";
    print "\n";
    print "OPTIONS\n";
    print "--units  or -u : specify the unit of energy. One of : mHa, meV, THz or cm-1.\n";
    print "--width  or -w : specify the width of the graph. defaults to 1. \n";
    print "--erange or -e : specify the range of the energy by the form [emin:emax].\n";
    print "--drange or -d : specify the range of the DOS by the form [dmin:dmax].\n";
    print "--title  or -t : specify the title of the graph. \n";
    print "--font   or -f : specify the fontsize to be used in the graph. defaults to 18.\n";
    print "--keep   or -k : specify this option in order to keep the intermediate data file.\n";
    print "--mono   or -m : specify this option in order to create monochrome figures.\n";
    print "--dinc         : specify the tic for the y-axis.\n";
    print "--einc         : specify the tic for the x-axis.\n";
}

sub parse_dos_data {
    open PHDOS,"<$dosfile";
    my @list=<PHDOS>;
    foreach my $line (@list) {
        if ($line =~ "#") {
            next;
        }
        $line = &strip($line);
        my @words=split(/\s+/,$line);
        if (scalar(@words)<8) {
            die("invalid line : $line read from file : $dosfile\n");
        }
        my @ed = ($words[1]*$unitfactor{$units},$words[4]/$unitfactor{$units});
        push(@dosdata,\@ed);
    }
    close(PHDOS);
}

sub write_dos_data {
    open PHDOS,">$tmpdosfile";
    foreach my $ddata (@dosdata) {
       print PHDOS "$ddata->[0] $ddata->[1]" . "\n"; 
    }
    close(PHDOS); 
}

sub run_gnuplot {
    my $gpi = "phonon_dos.gpi";
    open GP,">$gpi";
    if (defined($mono)) {
        print GP "set term post eps solid enhanced $fsize\n";
    } else {
        print GP "set term post eps color solid enhanced $fsize\n";
    }
    print GP "set output \"phonon_dos.eps\"\n";

    if($units eq 'meV') {
        print GP "set xlabel \"Frequency (meV)\"\n";
        print GP "set ylabel \"DOS (states/meV)\"\n";
    } elsif($units eq 'cm-1') {
        print GP "set xlabel \"Frequency (cm^{-1})\"\n";
        print GP "set ylabel \"DOS (states/cm^{-1})\"\n";
    } elsif($units eq 'THz') {
        print GP "set xlabel \"Frequency (THz)\"\n";
        print GP "set ylabel \"DOS (states/THz)\"\n";
    } else {
        print GP "set xlabel \"Frequency (mHa)\"\n";
        print GP "set ylabel \"DOS (states/mHa)\"\n";
    }

    if (defined($erange)) {
        print GP "set xrange $erange\n";
    }
    if (defined($drange)) {
        print GP "set yrange $drange\n";
    }
    if (defined($dinc)) {
        print GP "set ytic $dinc\n";
    }
    if (defined($einc)) {
        print GP "set xtic $einc\n";
    }
    if (defined($title)) {
        print GP "set title \"$title\"\n";
    }
    print GP "set size $width,1\n";
    print GP "plot \"$tmpdosfile\" w l not \n";
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
        `$rm $tmpdosfile\n`;
    }
}

