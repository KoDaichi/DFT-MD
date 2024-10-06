#!/usr/bin/perl -w

if(@ARGV < 1) {
   print "prep_strfrc.pl STRAIN INDEX_LIST\n";
   exit;
}

$strain = shift @ARGV;
print "STRAIN = $strain\n";
$index = shift @ARGV;
@index = split(' ',$index);
print "Index list = @index\n";

$num_strfrc = 0;
foreach $ie ( @index ) {
  foreach $pm ( 'p', 'm' ) {
    $e11="0.0"; $e22="0.0"; $e33="0.0";
    $e23="0.0"; $e32="0.0";
    $e31="0.0"; $e13="0.0";
    $e12="0.0"; $e21="0.0";
    if($pm eq 'p') {
      if($ie==1) {
        $e11 = $strain; 
      } elsif($ie==2) {
        $e22 = $strain; 
      } elsif($ie==3) {
        $e33 = $strain; 
      } elsif($ie==4) {
        $e23 = $strain/2; 
        $e32 = $strain/2; 
      } elsif($ie==5) {
        $e31 = $strain/2; 
        $e13 = $strain/2; 
      } elsif($ie==6) {
        $e12 = $strain/2; 
        $e21 = $strain/2; 
      }
    } elsif($pm eq 'm') {
      if($ie==1) {
        $e11 = -$strain; 
      } elsif($ie==2) {
        $e22 = -$strain; 
      } elsif($ie==3) {
        $e33 = -$strain; 
      } elsif($ie==4) {
        $e23 = -$strain/2; 
        $e32 = -$strain/2; 
      } elsif($ie==5) {
        $e31 = -$strain/2; 
        $e13 = -$strain/2; 
      } elsif($ie==6) {
        $e12 = -$strain/2; 
        $e21 = -$strain/2; 
      }
    }
    $scf_dir = "scf_${pm}${ie}";
    mkdir $scf_dir, 0777;
    `cp template_scf/file_names.data $scf_dir`;
    open(IN,"template_scf/nfinput.data");
    open(OUT,">$scf_dir/nfinput.data");
    while($line = <IN>) {
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
    }
    close(IN);
    close(OUT);
    $num_strfrc += 1;
  }
}

print "Number of strain-induced force calc. = $num_strfrc\n";
open(OUT,">exec_strfrc.pl");
print OUT <<HERE;
#!/usr/bin/perl -w

if(\@ARGV < 1) {
  print "exec_strfrc.pl PHASE PARALLEL\\n";
  print "Options:\\n";
  print " -arch={sr|primepower|vpp}\\n";
  print " -loadleveler\\n";
  exit;
}

\$mf = 'no';
if( -e 'machinefile' ) {
  \$mf = 'yes';
}

\$phase = shift \@ARGV;
\$num_proc = shift \@ARGV;

\$dummy=0;
\$machine='';
\$stop_parallel = 'no';
while ( \@ARGV > 0 ) {
  \$s = shift;
  if(\$s =~/-arch/) {
    (\$dummy,\$machine) = split('=',\$s);
  } elsif(\$s =~/-loadleveler/) {
    \$stop_parallel = 'yes';
  }
}

\@index = qw/ @index /;

\$id = 0;
\$num_id = \@index*2;

\$idd = 0;

loop:
# SCF
\$count = 0;
while( \$id < \$num_id ) {
\$count += 1;
\$id += 1;
\$idd += 1;

&set_id;

\$scf_dir = "scf_\${pm}\${ie}";

if(\$count==1) {
  open(OUT,">dirlist");
  print OUT "\$num_proc\\n";
  close(OUT);
  \$np = 0;
  \$nskip = 0;
}
if( -e \$scf_dir && -d \$scf_dir ) {
  if(\&check_scf_fin(\$scf_dir) eq 'no') {
    print " \$scf_dir calculating\\n";
    open(OUT,">>dirlist");
    print OUT "\$scf_dir\\n";
    close(OUT);
    \$np++;
  } else {
    print " \$scf_dir skipping\\n";
    \$nskip++;
  }
}
if( \$count == \$num_proc  || \$id == \$num_id ) {
  if(\$stop_parallel eq 'yes' && \$num_proc != \$np && \$nskip != \$num_proc) {
    print "This program stoped because you executed this program with -loadleveler option, and MPI parallelism is not the same to PARALLEL.\\n";
    print "Next MPI parallelism is \$np, then run this program with PARALLEL of \$num_proc and using \$np processes, and without the -loadleveler option.\\n";
    `rm dirlist`;
    exit;
  }
  if(\$np>0) {
    if(\$machine eq 'vpp') {
      `\$phase -np \$np`;
    } elsif(\$machine eq 'primepower') {
      `mpiexec -n \$np -mode limited \$phase`;
    } elsif(\$machine eq 'sr') {
      `mpirun -n \$np -np \$np \$phase`;
    } else {
      if(\$mf eq 'yes') {
        `mpirun -machinefile machinefile -np \$np \$phase`;
      } else {
        `mpirun -np \$np \$phase`;
      }
    }
  }
  `rm dirlist`;
  last;
}
#\$id -= \$count;
}
#\$id += \$count;
if( \$idd < \$num_id ) {
  goto loop;
}

opendir(DIR,"./");
\@files=readdir(DIR);
closedir(DIR);
open(OUT,">strfrc.data");
print OUT "$num_strfrc\\n";
foreach \$dir ( \@files ) {
  next if( !( -d \$dir ) || !((\$dir=~/^scf_p/) || (\$dir=~/^scf_m/)) );
  if( !(-f "\$dir/nfdynm.data") || (-f "\$dir/nfdynm.data" && -z "\$dir/nfdynm.data") ) {
    print "The calculation of \${dir} is not complelte.\n";
    exit;
  }
  if(\$dir =~/^scf_p/) {
    (\$dummy,\$ie) = split('p',\$dir);
    print OUT "\$ie $strain\\n";
  } elsif(\$dir =~/^scf_m/) {
    (\$dummy,\$ie) = split('m',\$dir);
    print OUT "\$ie -$strain\\n";
  }
  open(IN,"\$dir/nfdynm.data"); 
  while(\$line=<IN>) {
    if(\$line=~/ntyp/) {
      \@line = split('=',\$line);
      \$natom = \$line[2];
    } elsif(\$line=~/cps/) {
      for(\$i=0;\$i<\$natom;\$i++) {
        \@line = split(' ',<IN>);
        print OUT "\$line[4] \$line[5] \$line[6]\\n";
      } 
    }
  } 
  close(IN);
}
close(OUT);

# subroutine
sub set_id {
  \$i = 0;
  foreach \$iie ( \@index ) {
    foreach \$ipm ( 'p', 'm' ) {
      \$i += 1;
      if( \$id == \$i ) {
        \$ie = \$iie;
        \$pm = \$ipm;
        last;
      }
    }
  }
}

sub check_scf_fin(\$) {
  my (\$dir) = \@_;

  if( !(-f "\$dir/continue.data") ) {
    return 'no';
  }

  my \$iconv = 0;
  open(IN,"\$dir/continue.data");
  while(\$line=<IN>) {
    if(\$line =~/convergence/) {
      \$iconv = <IN>;
      chomp(\$iconv);
      last;
    }
  }
  close(IN);

  if(\$iconv == 0) {
    return 'no';
  } else {
    return 'yes';
  }
}

HERE
close(OUT);
chmod 0544,"exec_strfrc.pl";
