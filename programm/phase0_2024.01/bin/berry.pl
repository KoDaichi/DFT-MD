#!/usr/bin/perl -w

use strict;
use File::Copy 'copy';
use File::stat;
use threads;
use threads::shared;

sub trim($);

our %control_parameters : shared;
%control_parameters = (
"property"=>"zeff",
"cpumax"=>-1,
"stopcheck"=>10,
"template_scf"=>"template_scf",
"template_berry"=>"template_berry",
"mesh1"=>"",
"mesh2"=>"",
"mesh3"=>"",
"n1_1"=>4,
"n2_1"=>4,
"J_1"=>15,
"n1_2"=>4,
"n2_2"=>4,
"J_2"=>15,
"n1_3"=>4,
"n2_3"=>4,
"J_3"=>15,
"atom_list"=>undef,
"strain_list"=>undef,
"displacement"=>0.1,
"strain"=>0.01,
"scf_command"=>"mpirun phase",
"berry_command"=>"mpirun ekcal",
"ndir"=>1,
"np"=>1,
"ne"=>1,
"nk"=>1,
"ng"=>1,
"ng_b"=>1,
"ne_b"=>1,
"natm"=>undef,
"a_vector"=>undef,
"b_vector"=>undef,
"c_vector"=>undef,
"length_unit"=>'bohr',
"eval_zeff"=>1
);

our $clean;
our $mode="gendir";
our $control="control";
our $phonon_dir;
our $starttime=time;
our $in_termination : shared;
our @working_directories : shared;
our $term = undef;
our $bohr = 0.5291772480;

our $displacement_block="displacement{
    sw_displace_atom = on
    displaced_atom = <ATOM_ID>
    ux =  <Ux>
    uy =  <Uy>
    uz =  <Uz>
}
";

our $strain_block="strain{
    sw_strained_cell = ON
    e11 = <E11>, e22 = <E22>, e33 = <E33>
    e23 = <E23>, e32 = <E32>
    e31 = <E31>, e13 = <E13>
    e12 = <E12>, e21 = <E21>
}";

our $berry_block="berry_phase{
    sw_berry_phase = on
    g_index = <G_INDEX>
    mesh{ n1 = <MESH_N1>, n2 = <MESH_N2>, J = <MESH_J> }
}
";

our $zeff_ctrl_block="control{
    condition = initial
    max_iteration = 0
}
";

our $zeff_pp_block="postprocessing{
    polarization{
        sw_bp_property = on
        property = effective_charge
    }
}";

our $fchgt="F_CHGT   = '../<SCF_DIR>/nfchgt.data'";

our $focc="F_OCCMAT   = '../<SCF_DIR>/occmat.data'";

our $pawcntn="F_CNTN_BIN_PAW = '../<SCF_DIR>/continue_bin_paw.data'";

$in_termination = 0;

&main;

sub main {
    &parse_args;
    print "berry.pl : script to calculate the berry phase for the PHASE System\n";
    print "Copyright (c) Phase System Consortium \n";
    print 'E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp' . "\n";
    print "script start time : ".localtime(time)."\n";
    print "\n";
    print "-- parsing the control file --\n";
    &read_control;
    if($mode eq 'analyze') {
        &end;
    }
    print "-- generating directories --\n";
    (my $sdir, my $bdir) = &generate_directories;
    &dump_dirname;
    if($mode eq 'gendir'){
        &end("run this script with the '--mode=exec' option in order to actually perform the calculations.\n");
    }
    $term = threads->new(\&terminator);
    print "\n";
    print "-- doing SCF calculations --\n";
    &execute("target_directories_scf",$control_parameters{scf_command});
    if(!$in_termination) {
        print "-- done SCF calculations --\n";
        print "\n";
        if($control_parameters{property} eq 'strfrc') {
            &gen_strfrc_data($sdir);
        }
    }
    if(!$in_termination and 
       ($control_parameters{property} eq 'zeff' or 
        $control_parameters{property} eq 'piezo')) {
        print "-- doing berry calculations --\n";
        &execute("target_directories_berry",$control_parameters{berry_command});
        if(!$in_termination) {
            print "-- done berry calculations --\n";
            if(&cat_berry_data) {
                print "-- concatenated the obtained berry-phase data to : berry.data --\n";
                print "\n";
                print "done calculation of the berry-phase suited for $control_parameters{property}.\n";
            }
        }
    }
    if($control_parameters{property} eq 'zeff' and $control_parameters{eval_zeff}) {
        print "-- evaluating Zeff --\n";
#        &exec_zeff;
        &exec_zeff_via_dirlist;
        print "-- results --\n";
        &print_zeff(&get_newest("zeff","output[0-9][0-9][0-9]"));
    }
    &end;
}

sub end {
    my $narg = @_;
    if($narg>0) {
        print $_[0];
    }
    print "script end time : ".localtime(time)."\n";
    my $etime = time-$starttime;
    print "elapsed time : $etime (s)\n";
    if(defined($term)) {
        $term->detach;
    }
    exit;
}

sub clean {
    opendir DIR, ".";
    my @files = grep { !m/^(\.|\.\.)$/g } readdir DIR;
    closedir DIR;
    for my $file (@files) {
        if ( -d $file ) {
            if ( $file =~ "^berry_a+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^scf_a+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^scf_e+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^berry_e+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^scf_p+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^scf_m+" ) {
                 &rmrf($file);
            }
            if ( $file =~ "^zeff" ) {
                 &rmrf($file);
            }
        }
    }
    unlink "target_directories_berry";
    unlink "target_directories_scf";
    unlink "dirlist";
    unlink "berry.data";
}

sub rmrf {
    my $dir = $_[0];
    if(!-e $dir) {
        return;
    }
    opendir D, "$dir";
    my @files = grep { !m/^(\.|\.\.)$/g } readdir D;
    closedir D;
    for my $file (@files) {
        if( -d $file ){
            &rmrf($file);
        } else {
            unlink "$dir/$file";
        }
    }
    rmdir $dir;
}

sub gen_strfrc_data {
    my $num_strfrc = $_[0];
    my $strain = $control_parameters{strain};
    opendir(DIR,"./");
    my @files=readdir(DIR);
    closedir(DIR);
    open(OUT,">strfrc.data");
    print OUT "$num_strfrc\n";
    foreach my $dir ( @files ) {
        next if( !( -d $dir ) || !(($dir=~/^scf_p/) || ($dir=~/^scf_m/)) );
        if( !(-f "$dir/nfdynm.data") || (-f "$dir/nfdynm.data" && -z "$dir/nfdynm.data") ) {
            print "The calculation of ${dir} is not complelte.\n";
            exit;
        }
        if($dir =~/^scf_p/) {
            (my $dummy,my $ie) = split('p',$dir);
            print OUT "$ie $strain\n";
        } elsif($dir =~/^scf_m/) {
            (my $dummy,my $ie) = split('m',$dir);
            print OUT "$ie -$strain\n";
        }
        open(IN,"$dir/nfdynm.data"); 
        my $natom = 0;
        while(my $line=<IN>) {
            if($line=~/ntyp/) {
                my @words = split('=',$line);
                $natom = $words[2];
            } elsif($line=~/cps/) {
                for(my $i=0;$i<$natom;$i++) {
                    my @words = split(' ',<IN>);
                    print OUT "$words[4] $words[5] $words[6]\n";
                } 
            }
        } 
        close(IN);
    }
    close(OUT);
}

sub cat_berry_data {
    open(IN,"target_directories_berry") or &end("failed to open file : target_directories_berry\n");
    my @dirs = <IN>;
    my $ndir = @dirs;
    close(IN);
    open(OUT,">berry.data");
    print OUT "$ndir\n";
    close(OUT);
    foreach my $dir (@dirs) {
        chomp($dir);
        if (! -e "${dir}/berry.data") {
            print "ERROR! file does not exist : $dir/berry.data\n";
            unlink "berry.data";
            return 0;
        } else {
            open(IN,"$dir/berry.data");
            my @berry = <IN>;
            close(IN);
            open(OUT,">>berry.data");
            for my $b (@berry) {
                print OUT $b;
            }
            close(OUT);
        }
    }
    return 1;
}

sub terminator {
    for(;;){
        sleep $control_parameters{stopcheck};
        my $etime = time-$starttime;
        if ($control_parameters{cpumax}>0 and $etime>$control_parameters{cpumax}) {
            print "elapsed time ($etime s) exceeded the value specified by the [cpumax] variable ($control_parameters{cpumax} s)\n";
            print "terminating ASAP\n";
            $in_termination = 1;
        }
        if( -e 'stop' ) {
            print "detected the 'stop' file\n";
            print "terminating ASAP\n";
            $in_termination = 1;
            unlink 'stop';
        }
        if($in_termination){
            foreach my $dir (@working_directories) {
                open(ST,">$dir/nfstop.data");
                print ST "1\n";
                close(ST);
            }
            return;
        }
    }
}

sub parse_args {
    use Getopt::Std;
    use Getopt::Long;
    my $help;
    Getopt::Long::Configure('bundling');
    GetOptions(
    'm|mode=s'=>\$mode,
    'h|help'=> \$help,
    'c|clean'=>\$clean,
    'p|phonon_dir=s'=>\$phonon_dir
    );
    if(defined($clean)){
        &clean;
        exit;
    }
    if($mode ne 'gendir' and $mode ne 'exec' and $mode ne 'analyze'){
        &end("'mode' must be one of : 'gendir', 'exec' or 'analyze'");
    }
    my $len=@ARGV;
    if ($len==0 or defined($help)) {
        &print_usage;
        exit;
    }
    $control = $ARGV[0];
    if (! -e $control ){
        &end("control file : $control does not exist\n");
    }
}

sub check_parameters {
    if($control_parameters{property} ne "zeff"  and 
       $control_parameters{property} ne "piezo" and 
       $control_parameters{property} ne "strfrc") {
        &end("'property' must be one of : zeff, piezo or strfrc\n");
    }
    my $np = $control_parameters{np};
    my $nd = $control_parameters{ndir};
    my $ne = $control_parameters{ne};
    my $nk = $control_parameters{nk};
    my $ng = $control_parameters{ng};
    my $ng_b = $control_parameters{ng_b};
    my $ne_b = $control_parameters{ne_b};
    if( $np!=$nd*$ne*$nk*$ng ) {
        &end("np is not equal to nd x ne x nk x ng\n");
    }
    if( $np!=$nd*$ne_b*$ng_b ) {
        &end("np is not equal to nd x ne_b x ng_b\n");
    }

    if($control_parameters{property} eq "zeff" and !defined($control_parameters{atom_list})) {
        &end("atom_list must be defined when property = zeff\n");
    }
    if($control_parameters{property} eq "piezo" or $control_parameters{property} eq "strfrc") {
        if(!defined($control_parameters{strain_list})) {
            &end("strain_list must be defined when property = piezo or strfrc.\n");
        }
    }
}

sub estimate_mesh_parameters {
    if(!(defined($control_parameters{a_vector}) and 
         defined($control_parameters{b_vector}) and 
         defined($control_parameters{c_vector}))) {
        return;
    }
    my @avec = split(' ',trim($control_parameters{a_vector}));
    my @bvec = split(' ',trim($control_parameters{b_vector}));
    my @cvec = split(' ',trim($control_parameters{c_vector}));
    my $uni = lc($control_parameters{length_unit});
    @avec = (&convert($avec[0]),&convert($avec[1]),&convert($avec[2]));
    @bvec = (&convert($bvec[0]),&convert($bvec[1]),&convert($bvec[2]));
    @cvec = (&convert($cvec[0]),&convert($cvec[1]),&convert($cvec[2]));

    print "\n";
    print "unit cell (in bohr units)\n";
    print "a_vector : $avec[0] $avec[1] $avec[2] \n";
    print "b_vector : $bvec[0] $bvec[1] $bvec[2] \n";
    print "c_vector : $cvec[0] $cvec[1] $cvec[2] \n";

    my @a=();
    for(my $j=0 ; $j<3 ; $j++) {
        $a[0][$j] = $avec[$j];
        $a[1][$j] = $bvec[$j];
        $a[2][$j] = $cvec[$j];
    }

    my $deter = $a[0][0]*($a[1][1]*$a[2][2]-$a[2][1]*$a[1][2])
              + $a[1][0]*($a[2][1]*$a[0][2]-$a[0][1]*$a[2][2])
              + $a[2][0]*($a[0][1]*$a[1][2]-$a[1][1]*$a[0][2]);

    my @b=();
    $b[0][0]=($a[1][1]*$a[2][2]-$a[2][1]*$a[1][2])/$deter;
    $b[1][0]=($a[1][2]*$a[2][0]-$a[2][2]*$a[1][0])/$deter;
    $b[2][0]=($a[1][0]*$a[2][1]-$a[2][0]*$a[1][1])/$deter;
    $b[0][1]=($a[2][1]*$a[0][2]-$a[0][1]*$a[2][2])/$deter;
    $b[1][1]=($a[2][2]*$a[0][0]-$a[0][2]*$a[2][0])/$deter;
    $b[2][1]=($a[2][0]*$a[0][1]-$a[0][0]*$a[2][1])/$deter;
    $b[0][2]=($a[0][1]*$a[1][2]-$a[1][1]*$a[0][2])/$deter;
    $b[1][2]=($a[0][2]*$a[1][0]-$a[1][2]*$a[0][0])/$deter;
    $b[2][2]=($a[0][0]*$a[1][1]-$a[1][0]*$a[0][1])/$deter;

    my @b1;
    my @b2;
    my @b3;
    for(my $i=0 ; $i<3 ; $i++) {
        $b1[$i] = $b[$i][0];
        $b2[$i] = $b[$i][1];
        $b3[$i] = $b[$i][2];
    }
    my $c1;
    my $c2;
    my $dpara1;
    my $dpara2;
    my $dperp;
    my $n1;
    my $n2;
    my $J;
    my @mesh;
    ($c1, $c2)=&get_c1_c2(\@b2,\@b3,\@b1);
    $dpara1 = &norm($c1);
    $dpara2 = &norm($c2);
    $dperp  = &norm(\@b1);
    ($n1,$n2,$J) = &get_reference_mesh(($dpara1,$dpara2,$dperp));
    $mesh[0] = $n1;
    $mesh[1] = $n2;
    $mesh[2] = $J;
    print "\n";
    print "|b_para1|, |b_para2| and |b_perp| (in bohr^-1 units)\n";
    print "for reciprocal vector no. 1 : $dpara1, $dpara2, $dperp\n";

    ($c1, $c2)=&get_c1_c2(\@b1,\@b3,\@b2);
    $dpara1 = &norm($c1);
    $dpara2 = &norm($c2);
    $dperp  = &norm(\@b2);
    ($n1,$n2,$J) = &get_reference_mesh(($dpara1,$dpara2,$dperp));
    $mesh[3] = $n1;
    $mesh[4] = $n2;
    $mesh[5] = $J;
    print "for reciprocal vector no. 2 : $dpara1, $dpara2, $dperp\n";

    ($c1, $c2)=&get_c1_c2(\@b1,\@b2,\@b3);
    $dpara1 = &norm($c1);
    $dpara2 = &norm($c2);
    $dperp  = &norm(\@b3);
    ($n1,$n2,$J) = &get_reference_mesh(($dpara1,$dpara2,$dperp));
    $mesh[6] = $n1;
    $mesh[7] = $n2;
    $mesh[8] = $J;
    print "for reciprocal vector no. 3 : $dpara1, $dpara2, $dperp\n";
    print "\n";
    
    print "reference value for mesh parameters n1, n2 and J\n";
    print "for reciprocal vector no. 1 : $mesh[0], $mesh[1], $mesh[2]\n";
    print "for reciprocal vector no. 2 : $mesh[3], $mesh[4], $mesh[5]\n";
    print "for reciprocal vector no. 3 : $mesh[6], $mesh[7], $mesh[8]\n";
    print "\n";
}

sub get_reference_mesh {
    my $dpara1 = $_[0];
    my $dpara2 = $_[1];
    my $dpara3 = $_[2];
    my $n1 = int($dpara1/0.02);
    my $n2 = int($dpara2/0.02);
    my $J = int($dpara3/0.01);
    if($n1<1) {
        $n1 = 1;
    }
    if($n2<1) {
        $n2 = 1;
    }
    if($J<1) {
        $J = 1;
    }
    return ($n1,$n2,$J);
}

sub get_c1_c2 {
    my ($b1,$b2,$b3) = @_;
    my $bnorm2 = &norm($b3);
    $bnorm2 = $bnorm2*$bnorm2;
    my $dot1 = &dot_product($b1,$b3);
    my $dot2 = &dot_product($b2,$b3);

    my @c1 = (); 
    my @c2 = (); 
    for(my $i=0 ; $i<3 ; $i++) {
        $c1[$i] = @$b1[$i]-@$b3[$i]*$dot1/$bnorm2;
        $c2[$i] = @$b2[$i]-@$b3[$i]*$dot2/$bnorm2;
    }
    return (\@c1,\@c2);
}

sub dot_product {
    my ($vec_a, $vec_b) = @_;
    my $sum = 0;
    $sum += $vec_a->[$_] * $vec_b->[$_] for 0..$#$vec_a;
    return $sum;
}

sub norm {
    my $vec = $_[0];
    my $n=0.0;
    for(my $i=0 ; $i<3 ; $i++) {
        $n += @$vec[$i]*@$vec[$i];
    }
    return sqrt($n);
}

sub convert {
    my $org = $_[0];
    if($control_parameters{length_unit} eq 'bohr') {
        return $org;
    }
    if($control_parameters{length_unit} eq 'angstrom') {
        return $org/$bohr; 
    }
    if($control_parameters{length_unit} eq 'nm') {
        return 10*$org/$bohr;
    }
    &end("unknown length unit : $control_parameters{length_unit}\n");
}

sub read_control {
    open CTRL,"<$control" or 
      &end("failed to open control file : $control\n");
    while(my $line=<CTRL>){
        $line = trim($line);
        if($line =~ "^#"){
            next;
        }
        if(length($line)==0) {
            next;
        }
        my @words = split("#",$line);
        $line = $words[0];
        @words = split("=",$line,2);
        my $nwords = @words; 
        if ($nwords<2) {
            next;
        }
        my $key=trim($words[0]);
        if (exists $control_parameters{$key}){
            $control_parameters{$key}=trim($words[1]);
        }
    }
    close(CTRL);
    &check_parameters;
    &estimate_mesh_parameters;
    print "[parameters resolved from the control file]\n"; 
    print "property : $control_parameters{property}\n";
    print "cpumax   : $control_parameters{cpumax} (s)\n";
    print "check termination every  : $control_parameters{stopcheck} (s)\n";
    print "template scf directory   : $control_parameters{template_scf}\n"; print "template berry directory : $control_parameters{template_berry}\n"; 
    my @test1 = split(' ',$control_parameters{mesh1});
    my $n1 = @test1;
    if($n1 >= 3) {
        $control_parameters{n1_1} = $test1[0];
        $control_parameters{n2_1} = $test1[1];
        $control_parameters{J_1} = $test1[2];
    }
    my @test2 = split(' ',$control_parameters{mesh2});
    my $n2 = @test2;
    if($n2 >= 3) {
        $control_parameters{n1_2} = $test2[0];
        $control_parameters{n2_2} = $test2[1];
        $control_parameters{J_2} = $test2[2];
    }
    my @test3 = split(' ',$control_parameters{mesh3});
    my $n3 = @test3;
    if($n3 >= 3) {
        $control_parameters{n1_3} = $test3[0];
        $control_parameters{n2_3} = $test3[1];
        $control_parameters{J_3} = $test3[2];
    }
    if($control_parameters{property} eq "piezo" or $control_parameters{property} eq "zeff") {
        print "n1, n2 and J for the first  reciprocal lattice : $control_parameters{n1_1},$control_parameters{n2_1},$control_parameters{J_1}\n";
        print "n1, n2 and J for the second reciprocal lattice : $control_parameters{n1_2},$control_parameters{n2_2},$control_parameters{J_2}\n";
        print "n1, n2 and J for the third  reciprocal lattice : $control_parameters{n1_3},$control_parameters{n2_3},$control_parameters{J_3}\n";
    }
    if($control_parameters{property} eq "zeff") {
        print "atoms to be taken into account\n";
        print "$control_parameters{atom_list}\n";
        print "displacement : $control_parameters{displacement} (bohr)\n";
    } elsif ($control_parameters{property} eq "piezo" or $control_parameters{property} eq "strfrc") {
        print "components of the strain tensor to be taken into account\n";
        print "$control_parameters{strain_list}\n";
        print "strain : $control_parameters{strain}\n";
    }
    print "number of MPI processes : $control_parameters{np}\n";
    print "number of directories to be handled in parallel : $control_parameters{ndir}\n";
    print "k-point parallelization for the SCF calculation      : $control_parameters{nk}\n";
    print "band-parallelization for the SCF calculation         : $control_parameters{ne}\n";
    print "G-point parallelization for the SCF calculation      : $control_parameters{ng}\n";
    print "band-parallelization for the berry-phase calculation : $control_parameters{ne_b}\n";
    print "G-point parallelization for the berry-phase calculation : $control_parameters{ng_b}\n";
    print "SCF   command : $control_parameters{scf_command}\n";
    print "berry command : $control_parameters{berry_command}\n";
    print "(remark : NP will be replaced by np, NE will be replaced by ne, NK will be replaced by nk, NG_B will be replaced by ng_b and NE_B will be replaced by ne_b)\n";
    print "\n";
}

sub dump_dirname {
    if($control_parameters{property} eq 'zeff'){
        my $alist = $control_parameters{atom_list};
        my @atoms = split(' ',$alist);
        my $natm = @atoms;
        open(DIRLIST,">target_directories_scf");
        print DIRLIST "scf_a0\n";
        foreach my $aid (@atoms) {
            foreach my $uid ( 1, 2, 3 ) {
                print DIRLIST "scf_a${aid}_u${uid}\n";
            }
        }
        close(DIRLIST);
        open(DIRLIST,">target_directories_berry");
        print DIRLIST "berry_a0_g1\n";
        print DIRLIST "berry_a0_g2\n";
        print DIRLIST "berry_a0_g3\n";
        foreach my $aid (@atoms) {
            foreach my $uid ( 1, 2, 3 ) {
                foreach my $gid ( 1, 2, 3 ) {
                    print DIRLIST "berry_a${aid}_u${uid}_g${gid}\n";
                }
            }
        }
        close(DIRLIST);
    } elsif ($control_parameters{property} eq 'piezo' or $control_parameters{property} eq 'strfrc') {
        open(DIRLIST,">target_directories_scf");
        my $slist = $control_parameters{strain_list};
        my @s = split(' ',$slist);
        my $ns = @s;
        if($control_parameters{property} eq 'piezo') {
            print DIRLIST "scf_e0\n";
        }
        foreach my $sid (@s) {
            if($control_parameters{property} eq 'piezo'){
                print DIRLIST "scf_e${sid}\n";
            } elsif ($control_parameters{property} eq 'strfrc') {
                print DIRLIST "scf_p${sid}\n";
                print DIRLIST "scf_m${sid}\n";
            }
        }
        close(DIRLIST);
        if($control_parameters{property} eq 'piezo') {
            open(DIRLIST,">target_directories_berry");
            print DIRLIST "berry_e0_g1\n";
            print DIRLIST "berry_e0_g2\n";
            print DIRLIST "berry_e0_g3\n";
            foreach my $sid (@s) {
                foreach my $gid ( 1, 2, 3 ) {
                    print DIRLIST "berry_e${sid}_g${gid}\n";
                }
            }
            close(DIRLIST);
        }
    } 
}

sub generate_directories {
    my $fname = "file_names.data";
    my $scf = $control_parameters{template_scf};
    my $berry = $control_parameters{template_berry};
    my $sdir=0;
    my $bdir=0;
    if (! -e "$scf/$fname" ) {
        &end("$scf/$fname does not exist.\n");
    }
    if ($control_parameters{property} eq 'zeff' or 
        $control_parameters{property} eq 'piezo') {
        if (! -e "$berry/$fname" ) {
            &end("$berry/$fname does not exist.\n");
        }
    }
    my $alist = undef;
    my @atoms = undef;
    if($control_parameters{property} eq 'piezo' or 
       $control_parameters{property} eq 'strfrc') {
        $alist = $control_parameters{strain_list};
        @atoms = split(' ',$alist);
    } elsif ($control_parameters{property} eq 'zeff') {
        $alist = $control_parameters{atom_list};
        @atoms = split(' ',$alist);
    }
    my $natm = @atoms;
    if($natm==0) {
        if($control_parameters{property} eq 'zeff') {
            &end("atom list undefined\n");
        } elsif($control_parameters{property} eq 'piezo' or 
                $control_parameters{property} eq 'strfrc') {
            &end("strain list undefined\n");
        }
    }
    my @tmp = &generate_scf_directory(0);
    $sdir += $tmp[0];
    $bdir += $tmp[1];
    foreach my $id (@atoms) {
        @tmp = &generate_scf_directory($id);
        $sdir += $tmp[0];
        $bdir += $tmp[1];
    }
    print "number of SCF directories   : $sdir\n";
    print "number of berry directories : $bdir\n";
    return ($sdir,$bdir);
}

sub check_existense_of_tags {
    my $file = $_[0];
    open(IN,"$file") or &end("failed to open $file\n");
    my $target = "";
    while (my $line = <IN>) {
        my @words = split('!',$line);
        $target = $target . $words[0];
    }
    close(IN);
    my @check = $_[1];
    foreach my $ch (@check) {
        if (index($target,$ch)==-1) {
            &end("tag : $ch does not exist in file : $file\n");
        }
    }
}

sub generate_scf_directory {
    my $id = $_[0];
    my $tmpldir = $control_parameters{template_scf};
    my $scf_dir = "scf_a$id";
    if($control_parameters{property} eq 'piezo') {
        $scf_dir = "scf_e$id";
    }
    if($control_parameters{property} eq 'strfrc') {
        $scf_dir = "scf";
    }
    my $sdir = "";
    my $dircount=0;
    my $berrydircount=0;
    my $nfinp = &get_nfinp("$tmpldir/file_names.data");
    my @tags = ("<ATOM_ID>","<Ux>","<Uy>","<Uz>");
    if($control_parameters{property} eq 'piezo' or 
       $control_parameters{property} eq 'strfrc') {
        @tags = ("<E11>","<E22>","<E33>","<E23>","<E32>","<E31>","<E13>","<E21>","<E12>");
    }
#    &check_existense_of_tags("$tmpldir/$nfinp",@tags);
    my $e11="0.0"; my $e22="0.0"; my $e33="0.0";
    my $e23="0.0"; my $e32="0.0";
    my $e31="0.0"; my $e13="0.0";
    my $e12="0.0"; my $e21="0.0";
    my $strain = $control_parameters{strain};
    if (($control_parameters{property} eq 'zeff' and $id==0) or 
         $control_parameters{property} eq 'piezo') {
        $sdir = $scf_dir;
        if (! -e $sdir){
            mkdir $sdir;
        }
        if($control_parameters{property} eq 'piezo') {
            if($id==1) {
                $e11 = $strain; 
            } elsif($id==2) {
                $e22 = $strain; 
            } elsif($id==3) {
                $e33 = $strain; 
            } elsif($id==4) {
                $e23 = $strain/2; 
                $e32 = $strain/2; 
            } elsif($id==5) {
                $e31 = $strain/2; 
                $e13 = $strain/2; 
            } elsif($id==6) {
                $e12 = $strain/2; 
                $e21 = $strain/2; 
            }
        }
        $dircount++;
        copy "$tmpldir/file_names.data",$sdir or 
          &end("failed to copy $tmpldir/file_names.data to $sdir\n");
        open(OUT,">$sdir/$nfinp") or 
          &end("failed to open $sdir/$nfinp\n");

        my $line = get_input("$tmpldir/$nfinp");
        if($control_parameters{property} eq 'zeff'){
            $line = &strip_block($line,"displacement",1);
            $line = &insert_block($line,"atom_list",$displacement_block,2);
            $line =~s/\<ATOM_ID\>/$id/;
            $line =~s/\<Ux\>/0.0/;
            $line =~s/\<Uy\>/0.0/;
            $line =~s/\<Uz\>/0.0/;
        } elsif ($control_parameters{property} eq 'piezo') {
            $line = &strip_block($line,"strain",1);
            $line = &insert_block($line,"structure",$strain_block,1);
            $line =~s/\<E11\>/$e11/;
            $line =~s/\<E22\>/$e22/;
            $line =~s/\<E33\>/$e33/;
            $line =~s/\<E23\>/$e23/;
            $line =~s/\<E32\>/$e32/;
            $line =~s/\<E31\>/$e31/;
            $line =~s/\<E13\>/$e13/;
            $line =~s/\<E12\>/$e12/;
            $line =~s/\<E21\>/$e21/;
        }
        print OUT $line;
        close(OUT);
        if($control_parameters{property} eq 'zeff'){
            $berrydircount+=&generate_berry_directory($id);
        } elsif ($control_parameters{property} eq 'piezo') {
            $berrydircount+= &generate_berry_directory(
            $id,$e11,$e22,$e33,$e32,$e23,$e31,$e13,$e21,$e12);
        }
    } elsif ($control_parameters{property} eq 'zeff' and $id>0) {
        foreach my $iu ( 1, 2, 3 ) {
            $sdir = "${scf_dir}_u${iu}";
            if (! -e $sdir){
                mkdir $sdir;
            }
            $dircount++;
            copy "$tmpldir/file_names.data",$sdir or 
              &end("failed to copy $tmpldir/file_names.data to $sdir\n");
            my $line = get_input("$tmpldir/$nfinp");
            $line = &strip_block($line,"displacement",1);
            $line = &insert_block($line,"atom_list",$displacement_block,2);
            open(OUT,">$sdir/$nfinp") or 
              &end("failed to open $sdir/$nfinp\n");
            my $ux = 0.0; 
            my $uy = 0.0; 
            my $uz = 0.0; 
            if ($iu==1) {
                $ux = &convert($control_parameters{displacement});
            }
            if ($iu==2) {
                $uy = &convert($control_parameters{displacement});
            }
            if ($iu==3) {
                $uz = &convert($control_parameters{displacement});
            }
            $line =~s/\<ATOM_ID\>/$id/;
            $line =~s/\<Ux\>/$ux/;
            $line =~s/\<Uy\>/$uy/;
            $line =~s/\<Uz\>/$uz/;
            print OUT $line;
            close(OUT);
            $berrydircount+=&generate_berry_directory($id,$iu,$ux,$uy,$uz);
        }
    } elsif ($control_parameters{property} eq 'strfrc' and $id>0) {
        my @plusminus = (1.0, -1.0);
        foreach my $pm ( @plusminus ) {
            $sdir = "${scf_dir}_p${id}";
            if($pm<0) {
                $sdir = "${scf_dir}_m${id}";
            }
            if (! -e $sdir){
                mkdir $sdir;
            }
            $dircount++;
            copy "$tmpldir/file_names.data",$sdir or 
              &end("failed to copy $tmpldir/file_names.data to $sdir\n");
            my $line = get_input("$tmpldir/$nfinp");
            $line = &strip_block($line,"strain",1);
            $line = &insert_block($line,"structure",$strain_block,1);
            open(OUT,">$sdir/$nfinp") or 
              &end("failed to open $sdir/$nfinp\n");
            my $str = $pm*$strain;
            if ($id==1) {
                $e11 = $str; 
            } elsif($id==2) {
                $e22 = $str; 
            } elsif($id==3) {
                $e33 = $str; 
            } elsif($id==4) {
                $e23 = $str/2; 
                $e32 = $str/2; 
            } elsif($id==5) {
                $e31 = $str/2; 
                $e13 = $str/2; 
            } elsif($id==6) {
                $e12 = $str/2; 
                $e21 = $str/2; 
            }
            $line =~s/\<E11\>/$e11/;
            $line =~s/\<E22\>/$e22/;
            $line =~s/\<E33\>/$e33/;
            $line =~s/\<E23\>/$e23/;
            $line =~s/\<E32\>/$e32/;
            $line =~s/\<E31\>/$e31/;
            $line =~s/\<E13\>/$e13/;
            $line =~s/\<E12\>/$e12/;
            $line =~s/\<E21\>/$e21/;
            print OUT $line;
            close(OUT);
        }
    }
    return ($dircount,$berrydircount);
}

sub generate_berry_directory {
    my $id = $_[0];
    my $len = @_;
    my $uid;
    my $ux=0.0;
    my $uy=0.0;
    my $uz=0.0;
    my $dir=0;
    if($len==5){
        $uid = $_[1];
        $ux  = $_[2];
        $uy  = $_[3];
        $uz  = $_[4];
    }

    my $e11="0.0"; my $e22="0.0"; my $e33="0.0";
    my $e23="0.0"; my $e32="0.0";
    my $e31="0.0"; my $e13="0.0";
    my $e12="0.0"; my $e21="0.0";
    if($len==10){
        $e11 = $_[1];
        $e22 = $_[2];
        $e33 = $_[3];
        $e32 = $_[4];
        $e23 = $_[5];
        $e31 = $_[6];
        $e13 = $_[7];
        $e21 = $_[8];
        $e12 = $_[9];
    }

    my $scf_dir = "scf_a${id}";
    if(defined($uid)) {
        $scf_dir = "${scf_dir}_u${uid}";
    }
    if($control_parameters{property} eq 'piezo') {
        $scf_dir = "scf_e${id}";
    }
    my $tmpldir = $control_parameters{template_berry};
    my $berry_dir = "berry_a${id}";
    if(defined($uid)) {
        $berry_dir = "${berry_dir}_u${uid}";
    }
    if($control_parameters{property} eq 'piezo') {
        $berry_dir = "berry_e${id}";
    }
    my $nfinp = &get_nfinp("$tmpldir/file_names.data");
    my @tags = ("<ATOM_ID>","<Ux>","<Uy>","<Uz>","<G_INDEX>","<MESH_N1>","<MESH_N2>","<MESH_J>");
    if($control_parameters{property} eq 'piezo') {
        @tags = ("<E11>","<E22>","<E33>","<E23>","<E32>","<E31>","<E13>","<E21>","<E12>");
    }
    foreach my $ig ( 1, 2, 3 ) {
        my $bdir = "${berry_dir}_g${ig}";
        if (! -e $bdir) {
            mkdir $bdir;
        }
        $dir += 1;

        open(IN,"$tmpldir/file_names.data") or 
          &end("failed to open $tmpldir/file_names.data\n");

        open(OUT,">$bdir/file_names.data") or 
          &end("failed to open $bdir/file_names.data\n");
        while(my $line = <IN>) {
            print OUT $line;
        }
        close(IN);
        close(OUT);

        my $fdata = &edit_file_names_data("$bdir/file_names.data");
        my @fdata = split("\n",$fdata);
        open(OUT,">$bdir/file_names.data") or 
          &end("failed to open $bdir/file_names.data\n");
        foreach my $line (@fdata) {
            $line =~s/\<SCF_DIR\>/$scf_dir/;
            print OUT $line . "\n";
        }
        close(OUT);

        open(OUT,">$bdir/$nfinp") or 
          &end("failed to open $bdir/$nfinp\n");
        my $n1 = $control_parameters{n1_1};
        my $n2 = $control_parameters{n2_1};
        my $J  = $control_parameters{J_1};
        if($ig==2) {
           $n1 = $control_parameters{n1_2};
           $n2 = $control_parameters{n2_2};
           $J  = $control_parameters{J_2};
        }
        if($ig==3) {
           $n1 = $control_parameters{n1_3};
           $n2 = $control_parameters{n2_3};
           $J  = $control_parameters{J_3};
        }
        my $line = get_input("$tmpldir/$nfinp");
        if($control_parameters{property} eq 'zeff') {
            $line = &strip_block($line,"displacement",1);
            $line = &insert_block($line,"atom_list",$displacement_block,2);
        } elsif ($control_parameters{property} eq 'piezo') {
            $line = &strip_block($line,"strain",1);
            $line = &insert_block($line,"structure",$strain_block,1);
        }
        $line = &strip_block($line,"berry_phase",0);
        $line = &insert_block($line,"",$berry_block,0);
        if($control_parameters{property} eq 'zeff') {
            $line =~s/\<ATOM_ID\>/$id/;
            $line =~s/\<Ux\>/$ux/;
            $line =~s/\<Uy\>/$uy/;
            $line =~s/\<Uz\>/$uz/;
            $line =~s/\<G_INDEX\>/$ig/;
            $line =~s/\<MESH_N1\>/$n1/;
            $line =~s/\<MESH_N2\>/$n2/;
            $line =~s/\<MESH_J\>/$J/;
        } elsif ($control_parameters{property} eq 'piezo') {
            $line =~s/\<E11\>/$e11/;
            $line =~s/\<E22\>/$e22/;
            $line =~s/\<E33\>/$e33/;
            $line =~s/\<E23\>/$e23/;
            $line =~s/\<E32\>/$e32/;
            $line =~s/\<E31\>/$e31/;
            $line =~s/\<E13\>/$e13/;
            $line =~s/\<E12\>/$e12/;
            $line =~s/\<E21\>/$e21/;
            $line =~s/\<G_INDEX\>/$ig/;
            $line =~s/\<MESH_N1\>/$n1/;
            $line =~s/\<MESH_N2\>/$n2/;
            $line =~s/\<MESH_J\>/$J/;
        }
        print OUT $line;
        close(OUT);
    }
    return $dir;
}

sub get_nfinp {
    my $fname = $_[0];
    open(IN,"$fname") or &end("failed to open $fname\n");
    while(my $line=<IN>) {
        $line = trim($line);
        if ($line =~ "F_INP") {
            my @words = split("=",$line);
            my $len = @words;
            if($len==2) {
                my $ret = trim($words[1]);
                my $l = length($ret);
                if ($l>2) {
                    close(IN);
                    return substr($ret,1,-1);
                }
            }
        }
    }
    close(IN);
    return "nfinp.data";
}

sub execute {
    my $target_dir = $_[0];
    my $command    = $_[1];
    open(IN,$target_dir) or &end("failed to open file : $target_dir\n");
    my @dirs = <IN>;
    close(IN);
    my $nparallel = $control_parameters{ndir};
    if ($nparallel<1) {
        &end("variable ndir must be larger than or equal to 1\n");
    }
    my $ndir = @dirs;
    if($ndir<1) {
        &end("number of directories registered in $target_dir is less than 1\n");
    }
    my $ncal = int($ndir/$nparallel);
    my $ndir_per_cal=0;
    if ($ncal>0) {
        $ndir_per_cal = $nparallel;
    }
    my $amari = $ndir % $nparallel;
    for (my $i=0 ; $i<$ncal ; $i++) {
        my $comm = $command;
        $comm =~s/NP/$control_parameters{np}/;
        $comm =~s/NE_B/$control_parameters{ne_b}/;
        $comm =~s/NG_B/$control_parameters{ng_b}/;
        $comm =~s/NE/$control_parameters{ne}/;
        $comm =~s/NK/$control_parameters{nk}/;
        $comm =~s/NG/$control_parameters{ng}/;
        if ($ndir_per_cal==1) {
            print "running [$comm] under the following directories\n";
        } else {
            print "running [$comm] under the following directories\n";
        }
        open(DIRLIST,">dirlist") or 
          &end("failed to create the dirlist file\n"); 
        print DIRLIST "$nparallel\n";
        my $geta = $i*$nparallel;
        @working_directories=();
        my $converged_directories = 0;
        for (my $j=0 ; $j<$ndir_per_cal ; $j++){
            print DIRLIST "$dirs[$geta+$j]";
            print $dirs[$geta+$j];
            chomp($dirs[$geta+$j]);
            $working_directories[$j] = $dirs[$geta+$j];
            if(&is_converged($dirs[$geta+$j])) {
                $converged_directories++;
            }
        }
        close DIRLIST;
        if($converged_directories!=$ndir_per_cal) {
            &run($comm);
            if($in_termination) {
                return;
            }
        }
    }
    if($amari>0) {
        for (my $i=0 ; $i<$amari ; $i++) {
            my $comm = $command;
            $comm =~s/NP/$control_parameters{np}/;
            $comm =~s/NE_B//g;
            $comm =~s/NE//g;
            $comm =~s/NK//g;
            $comm =~s/ne//g;
            $comm =~s/nk//g;
            $comm =~s/ng//g;
            $comm =~s/=//g;
            print "running [$comm] under the following directory\n";
            open(DIRLIST,">dirlist") or 
              &end("failed to create the dirlist file\n"); 
            @working_directories=();
            print DIRLIST "1\n";
            print DIRLIST "$dirs[$ncal*$nparallel+$i]";
            print $dirs[$ncal*$nparallel+$i];
            chomp($dirs[$ncal*$nparallel+$i]);
            $working_directories[0] = $dirs[$ncal*$nparallel+$i];
            close(DIRLIST);
            if(!&is_converged($dirs[$ncal*$nparallel+$i])){
                &run($comm);
                if($in_termination) {
                    return;
                }
            }
        }
    }
    unlink "dirlist";
}

sub is_converged {
    my $dir = $_[0];
    if ( !-e "$dir/continue.data" ){
        return 0;
    }
    open(IN,"$dir/continue.data") or 
        &end("failed to open : $dir/continue.data");
    while(my $line=<IN>) {
        chomp($line);
        if ($line=~"convergence") {
            $line = <IN>;
            $line = trim($line);
            chomp($line);
            if ($line=~"2") {
                print "calculation under $dir seems to be converged.\n";
                close(IN);
                return 1;
            }
        }
    }
    close(IN);
    return 0;
}

sub run {
    my $command = $_[0];
    my $start = time;
    my $ret = system($command);
    if( $ret!=0 ) {
        &end("encountered error while executing : [$command]\n");
    }
    my $etime = time-$start;
    my $eltime = time-$starttime;
    print "time spent in this calculation : $etime (s), total time : $eltime (s)\n";
}

sub print_usage {
    print "Usage : berry.pl control [OPTIONS]\n";
}

sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub get_ws {
    my $nws = $_[0];
    my $ws = "";
    for(my $i=0 ; $i<$nws ; $i++) {
        $ws = $ws . "    ";
    }
    return $ws;
}

sub strip_block {
    my @lines = split("\n",$_[0]);
    my $tag = $_[1] . "{";
    my $print_ws = $_[2];
    my $count = 0;
    my $depth = 0;
    my $target_depth = -1;
    my $bl = "";
    my $ret = "";
    my $in_bl=0;
    my $in_target_block = 0;
    foreach my $line ( @lines ) {
        if( $line =~ ('\w*{\w*')) {
            if( lc($line) !~ ('\A' . quotemeta($tag)) && !$in_target_block) {
                if($print_ws) {
                    $ret = $ret . &get_ws($depth) . $line . "\n";
                } else {
                    $ret = $ret . $line . "\n";
                }
            }
            $depth += 1;
            $in_bl= 1;
        }
        if( $line =~ ('\w*}\w*')) {
            if($target_depth != $depth && !$in_target_block) {
                if($print_ws) {
                    $ret = $ret . &get_ws($depth-1) . $line . "\n";
                } else {
                    $ret = $ret . $line . "\n";
                }
            }
            if($target_depth == $depth) {
                $bl = $bl . "}\n";
                $target_depth = -1;
                $in_target_block = 0;
            }
            $depth -= 1;
            $in_bl= 1;
        }
        if( lc($line) =~ ('\A' . quotemeta($tag))) {
            $target_depth = $depth;
            $in_target_block = 1;
        }
        if($target_depth >0) {
            $bl = $bl . $line . "\n";
        } else {
          if(!$in_bl) {
              if($print_ws) {
                  $ret = $ret . &get_ws($depth) . $line . "\n";
              } else {
                  $ret = $ret . $line . "\n";
              }
          }
        }
        $in_bl= 0;
    }

    return $ret;
}

sub insert_block {
    my @lines = split("\n",$_[0]);
    my $tag = $_[1] . "{";
    my $blk = $_[2];
    my $target_depth = $_[3];
    my $depth = 0;
    my $ret = "";
    if( $target_depth == 0 ) {
        $ret = $_[0] . $_[2];
        return $ret;
    }
    foreach my $line ( @lines ) {
        if( $line =~ ('\w*{\w*')) {
            $depth += 1;
        }
        if( $line =~ ('\w*}\w*')) {
            $depth -= 1;
        }
        $ret = $ret . $line . "\n";
        if( trim(lc($line)) =~ ('\A' . quotemeta($tag)) && $depth == $target_depth) {
            my @blks = split("\n",$blk);
            foreach my $b (@blks) {
                $ret = $ret . get_ws($depth) . $b . "\n";
            }
        }
    }
    return $ret;
}

sub get_input {
    my $file = $_[0];
    open(IN,"$file") or &end("failed to open $file\n");
    my $target = "";
    while (my $line = <IN>) {
        my $word = trim($line);
        if(index($word,'!#') == -1) {
            my @words = split('!',$word);
            if($#words>0) {
                $word = $words[0];
            }
        }
        my $poss = index($word,"{");
        my $pose = index($word,"}");
        my @wordss = split("{",$word);
        if($poss>0) {
            $wordss[0] = $wordss[0] . "{";
        }
        my $theline = "";
        foreach my $w (@wordss) {
            $theline = $theline . $w . "\n";
        }
        my @wordse = split("}",$theline);
        if($pose>0) {
           $theline = "";
           foreach my $w (@wordse) {
               $theline = $theline . $w;
           }
           $theline = $theline . "}\n";
        }

        $target = $target . $theline;
    }
    close(IN);
    return $target;
}

sub edit_file_names_data {
   my $fname = $_[0];
   open(IN,"$fname");
   my $acc = "";
   while (my $line=<IN>) {
     $line = trim($line);
     if (trim(lc($line)) !~ ('\A' . "f_chgt")   &&
         trim(lc($line)) !~ ('\A' . "f_occmat") &&
         trim(lc($line)) !~ ('\A' . "f_cntn_bin_paw")) {
         $acc = $acc . $line . "\n"; 
     }
   }
   close(IN);
   my @lines = split("\n",$acc);
   $acc = "";
   foreach my $line (@lines) {
       if( index($line,'fnames') != -1 ) {
           $acc = $acc . "&fnames" . "\n" . $fchgt . "\n";
           $acc = $acc . $focc . "\n";
           $acc = $acc . $pawcntn . "\n";
       } else {
           $acc = $acc . $line . "\n";
       }
   }
   return $acc;
}

sub exec_zeff_via_dirlist {
    my $zeff = "zeff";
    my $tmpldir = "scf_a0";
    my $berrydata = "berry.data";
    if(! -e $zeff) {
        mkdir $zeff;
    }
    my $nfinp = &get_nfinp("$tmpldir/file_names.data");
    copy "$tmpldir/file_names.data",$zeff or &end("failed to copy $tmpldir/file_names.data");
    copy "$tmpldir/$nfinp",$zeff or &end("failed to copy $tmpldir/file_names.data");
    copy "$berrydata",$zeff or &end("failed to copy $berrydata");
    my $line = get_input("$tmpldir/$nfinp");
    $line = &strip_block($line,"displacement",1);
    $line = &strip_block($line,"postprocessing",0);
    $line = &strip_block($line,"control",0);
    $line = &insert_block($line,"",$zeff_ctrl_block,0);
    $line = &insert_block($line,"",$zeff_pp_block,0);
    $line = $line . "\n";
    open(OUT,">$zeff/$nfinp");
    print OUT $line;
    close(OUT);

    open(DIRLIST,">dirlist") or 
      &end("failed to create the dirlist file\n"); 
    print DIRLIST "1\n";
    print DIRLIST "$zeff";
    close(DIRLIST);

    my $phase = $control_parameters{scf_command};
    $phase =~s/NP/1/;
    $phase =~s/NE/1/;
    $phase =~s/NK/1/;
    $phase =~s/NG/1/;
    &run($phase);
    unlink "dirlist";
}

sub exec_zeff {
    my $zeff = "zeff";
    my $tmpldir = "scf_a0";
    my $berrydata = "berry.data";
    if(! -e $zeff) {
        mkdir $zeff;
    }
    my $nfinp = &get_nfinp("$tmpldir/file_names.data");
    copy "$tmpldir/file_names.data",$zeff or &end("failed to copy $tmpldir/file_names.data");
    copy "$tmpldir/$nfinp",$zeff or &end("failed to copy $tmpldir/file_names.data");
    copy "$berrydata",$zeff or &end("failed to copy $berrydata");
    my $line = get_input("$tmpldir/$nfinp");
    $line = &strip_block($line,"displacement",1);
    $line = &strip_block($line,"postprocessing",0);
    $line = &strip_block($line,"control",0);
    $line = &insert_block($line,"",$zeff_ctrl_block,0);
    $line = &insert_block($line,"",$zeff_pp_block,0);
    $line = $line . "\n";
    open(OUT,">$zeff/$nfinp");
    print OUT $line;
    close(OUT);

    chdir($zeff);
    my $phase = $control_parameters{scf_command};
    $phase =~s/NP/1/;
    $phase =~s/NE/1/;
    $phase =~s/NK/1/;
    $phase =~s/NG/1/;
    &run($phase);
    chdir("..");
}

sub get_newest {
    my $dir = $_[0];
    my $fil = $_[1];
    opendir DIR, "$dir";
    my $newest = "";
    my $max = -1;
    while (defined (my $file = readdir(DIR))) {
        if($file ne "." && $file ne ".." && $file =~ $fil) {
            if(stat("$dir/$file")->mtime>$max) {
                $newest = $file;
                $max = stat("$dir/$file")->mtime;
            }
        }
    }
    return "$dir/$newest";
}

sub print_zeff {
    my $start_tag = "--- Calculated effective charges ---";
    my $end_tag = "F_EFFCHG";
    my $file = $_[0];
    open(IN,$file);
    my $in_berry = 0;
    while(my $line=<IN>) {
      if(trim($line) =~ ('\A' . $start_tag)) {
          $in_berry = 1;
      }
      if(trim($line) =~ ('\A' . $end_tag)) {
          $in_berry = 0;
      }
      if(!$in_berry) {
          next;
      }
      print $line;
    }
    close(IN);
}

#sub get_description_from_dirname {
#    my $dirname = $_[0];
#    if ( $dirname =~ ('\A' . "berry_")) {
#        my @a0 = split("berry_",$dirname);
#        my @atom = split(@a0);
#    }
#}

