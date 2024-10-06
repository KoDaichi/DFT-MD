#!/usr/bin/perl -w

if(@ARGV < 1) {
   print "prep_piezo.pl STRAIN INDEX_LIST MESH1 MESH2 MESH3\n";
   exit;
}

$strain = shift @ARGV;
print "STRAIN = $strain\n";
$index = shift @ARGV;
@index = split(' ',$index);
print "Index list = @index\n";
$mesh1 = shift @ARGV;
@mesh1 = split(' ',$mesh1);
print "Mesh_1 = @mesh1\n";
$mesh2 = shift @ARGV;
@mesh2 = split(' ',$mesh2);
print "Mesh_2 = @mesh2\n";
$mesh3 = shift @ARGV;
@mesh3 = split(' ',$mesh3);
print "Mesh_3 = @mesh3\n";

unshift @index, 0;

$num_berry = 0;
foreach $ie ( @index ) {
  $e11="0.0"; $e22="0.0"; $e33="0.0";
  $e23="0.0"; $e32="0.0";
  $e31="0.0"; $e13="0.0";
  $e12="0.0"; $e21="0.0";
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
  $scf_dir = "scf_e${ie}";
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
  foreach $ig ( 1, 2, 3 ) {
    if( $ig == 1 ) {
      $n1 = $mesh1[0]; $n2 = $mesh1[1]; $j = $mesh1[2];
    } elsif( $ig == 2 ) {
      $n1 = $mesh2[0]; $n2 = $mesh2[1]; $j = $mesh2[2];
    } elsif( $ig == 3 ) {
      $n1 = $mesh3[0]; $n2 = $mesh3[1]; $j = $mesh3[2];
    }
    $berry_dir = "berry_e${ie}_g${ig}";
    mkdir $berry_dir, 0777;
    open(IN,"template_berry/file_names.data");
    open(OUT,">$berry_dir/file_names.data");
    while($line = <IN>) {
      $line =~s/\<SCF_DIR\>/$scf_dir/;
      print OUT $line;
    }
    close(IN);
    close(OUT);
    open(IN,"template_berry/nfinput.data");
    open(OUT,">$berry_dir/nfinput.data");
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
      $line =~s/\<G_INDEX\>/$ig/;
      $line =~s/\<MESH_N1\>/$n1/;
      $line =~s/\<MESH_N2\>/$n2/;
      $line =~s/\<MESH_J\>/$j/;
      print OUT $line;
    }
    close(IN);
    close(OUT);
    $num_berry += 1;
  }
}

print "Number of berry calc. = $num_berry\n";
open(OUT,">exec_piezo.pl");
print OUT <<HERE;
#!/usr/bin/perl -w

if(\@ARGV < 1) {
  print "exec_piezo.pl PHASE EKCAL PARALLEL\\n";
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
\$ekcal = shift \@ARGV;
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
\$num_id = \@index;

\$idd = 0;

loop:
# SCF
\$count = 0;
while( \$id < \$num_id ) {
\$count += 1;
\$id += 1;
\$idd += 1;

&set_id;

\$scf_dir = "scf_e\${ie}";

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
}
\$id -= \$count;

# Berry
foreach \$ig ( 1, 2, 3 ) {
\$count = 0;
while( \$id < \$num_id ) {
\$count += 1;
\$id += 1;

&set_id;

\$berry_dir = "berry_e\${ie}_g\${ig}";
if(\$count==1) {
  open(OUT,">dirlist");
  print OUT "\$num_proc\\n";
  close(OUT);
  \$np = 0;
}
if( -e \$berry_dir && -d \$berry_dir ) {
  if(\&check_berry_fin(\$berry_dir) eq 'no') {
    print " \$berry_dir calculating\\n";
    open(OUT,">>dirlist");
    print OUT "\$berry_dir\\n";
    close(OUT);
    \$np++;
  } else {
    print " \$berry_dir skipping\\n";
  }
}
if( \$count == \$num_proc  || \$id == \$num_id ) {
  if(\$np>0) {
    if(\$machine eq 'vpp') {
      `\$ekcal -np \$np`;
    } elsif(\$machine eq 'primepower') {
      `mpiexec -n \$np -mode limited \$ekcal`;
    } elsif(\$machine eq 'sr') {
      `mpirun -n \$np -np \$np \$ekcal`;
    } else {
      if(\$mf eq 'yes') {
        `mpirun -machinefile machinefile -np \$np \$ekcal`;
      } else {
        `mpirun -np \$np \$ekcal`;
      }
    }
  }
  `rm dirlist`;
  last;
}
}
\$id -= \$count;
}
\$id += \$count;
if( \$idd < \$num_id ) {
  goto loop;
}

opendir(DIR,"./");
\@files=readdir(DIR);
closedir(DIR);
open(OUT,">berry.data");
print OUT "$num_berry\\n";
foreach \$dir ( \@files ) {
  next if( !( -d \$dir ) || !(\$dir=~/^berry_e/) );
  if( !(-f "\$dir/berry.data") || (-f "\$dir/berry.data" && -z "\$dir/berry.data") ) {
    print "The calculation of \${dir} is not complelte.\n";
    exit;
  }
  open(IN,"\$dir/berry.data"); 
  print OUT <IN>;
  close(IN);
}
close(OUT);

# subroutine
sub set_id {
  \$i = 0;
  foreach \$iie ( \@index ) {
    \$i += 1;
    if( \$id == \$i ) {
      \$ie = \$iie;
      last;
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

sub check_berry_fin(\$) {
  my (\$dir) = \@_;

  if( !( -f "\$dir/berry.data" ) ) {
    return 'no';
  } elsif( -z "\$dir/berry.data" ) {
    return 'no';
  } else {
    return 'yes';
  }
}
HERE
close(OUT);
chmod 0544,"exec_piezo.pl";
