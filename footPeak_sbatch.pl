#!/usr/bin/perl
# footLoop pipeline Version 1.2
# Copyright (C) 2019 Stella Regina Hartono
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The license can be found at https://www.gnu.org/licenses/gpl-3.0.en.html. 
# By downloading or using this software, you agree to the terms and conditions of the license. 

use strict; use warnings; use Getopt::Std; use Time::HiRes; use Benchmark qw(:all); use Benchmark ':hireswallclock'; use Carp; use Thread; use Thread::Queue; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_g $opt_p $opt_d $opt_s $opt_k $opt_K $opt_n $opt_h $opt_t $opt_w $opt_L $opt_l $opt_o $opt_A $opt_G $opt_Z $opt_O $opt_0 $opt_J $opt_F);
getopts("vg:p:d:s:k:K:n:ht:w:l:A:o:G:L:ZO:0J:F");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

#   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
  # $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
#$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
}

use myFootLib;
use FAlite;
use footPeakAddon;
use footLoop_sbatch;

my $homedir = $ENV{"HOME"};

my ($footLoop_script_folder, $version, $md5script) = check_software();
my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
$version_small = $version if not defined $version_small;

##################
# 0. Check Sanity #
##################

my $date = getDate();
my $uuid = getuuid();
my $numThread = 1;
my ($footLoopFolder) = $opt_n;
my ($usage) = check_sanity($footLoopFolder);
my $logFile = "$footLoopFolder/logFile.txt";
print "\nfootLooplogFile:\n$LCY$logFile$N\n\n";
my ($opts, $outLog) = parse_footLoop_logFile($logFile, $date, $uuid, $footLoopFolder, $version_small);
$opts->{label}    = defined $opt_L ? $opt_L : $opts->{label};
my $label         = $opts->{label};
my $seqFile       = $opts->{seqFile};
my $indexFile		= $opts->{geneIndexFile};
my $origFolder 	= $opts->{origDir};
my $minLen			= $opts->{l};
my $x 				= $opts->{x};
my $y 				= $opts->{y};
my $threshold 		= $opts->{t}; $threshold /= 100 if $threshold > 1;
my $window 			= $opts->{w};
my $min 				= $opts->{d};
my $groupsize 		= $opts->{s};
my $minDis 			= $opts->{k};
my $NDir       	= $opts->{n};
my $outDir        = $opts->{o};
my $genewant      = $opts->{G}; $genewant = -1 if $genewant eq "FALSE";

#########################
# 1. Define Input/Names #
#########################

$label = getFilename($opts->{BAMFile}, "full") if not defined $label;
if (not -e "$outDir/.LABEL") {
	open (my $outLabel, ">", "$outDir/.LABEL") or die "Failed to write to $outDir/.LABEL: $!\n";
	print $outLabel $label;
	close $outLabel;
}
system("/bin/cp $seqFile $outDir");
LOG($outLog, "seqFile=$seqFile\n");

makedir("$outDir/.CALL") if not -d "$outDir/.CALL";

LOG($outLog, $usage . "faFile=$seqFile\nindexFile=$indexFile\norigFolder=$origFolder\nNDir=$NDir\noutDir=$outDir\n") and die unless defined $seqFile and defined $indexFile and defined $origFolder and -e $seqFile and -e $indexFile and -e $origFolder and defined $NDir and -d $NDir and defined $outDir;

if (not -d $outDir) {
	system("mkdir -p $outDir") or LOG($outLog, "Cannot create directory (-n) $outDir: $!\n\n") and die;
}

my $outsbatchDir = "$outDir/.footPeak_sbatch/";
system("mkdir -p $outDir/.footPeak_sbatch/");

# Get Real Coordinate of each Gene
my $SEQ = parse_indexFile_and_seqFile($indexFile, $seqFile, $outLog);

# Get origFiles from footLoop folder
my @origFile = <$origFolder/origFiles/*.filtered.gz>;
LOG($outLog, "Found $LGN" . scalar(@origFile) . "$N *.filtered.gz files!\n");

open (my $outLogAddon, ">", "$outDir/footLoop_addition_logFile.txt") or DIELOG($outLog, "Failed to write to $outDir/footLoop_addition_logFile.txt: $!\n");
close $outLogAddon;
my $out;
my $R;
my @total_Rscript = <$origFolder/*.R>;
my @total_png = <$origFolder/*.png>;
my %genes;

for (my $i = 0; $i < @origFile; $i++) {
	my ($peakFolder, $peakFilename) = getFilename($origFile[$i], "folderfull");
	$peakFilename =~ s/.filtered.gz$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	$genes{uc($gene)} ++;
	LOG($outLog, "Example: $LCY$origFile[$i]\t$gene\t$strand$N\n\n") if $i eq 0;
}

my $totalorigFile = @origFile;
foreach my $origFile (@origFile) {
}
LOG($outLog, "\n\n--------------$YW\n1. footPeak.pl processing $LGN" . scalar(keys %genes) . "$N genes (total file = $LGN" . scalar(@origFile) . "$N)\n\n");

my ($max_parallel_run) = defined $opt_J ? $opt_J : 1;
#my $cmd = "$footLoop_script_folder/footPeak_sbatch_2.pl $indexFile $seqFile FNINDICE FILENAME $outDir $window $threshold $totalorigFile $label $genewant $minDis $minLen $version_small";
my $cmd = "$footLoop_script_folder/footPeak_sbatch_3.pl $indexFile $seqFile FNINDICE FILENAME $outDir $window $threshold $totalorigFile $label $genewant $minDis $minLen $version_small";
my $force_sbatch = $opt_F;

#my @origFilesmall = @origFile[0..1];
makedir("$outDir/PEAKS_GENOME/") if not -d "$outDir/PEAKS_GENOME/";
makedir("$outDir/PEAKS_LOCAL/") if not -d "$outDir/PEAKS_LOCAL/";
footLoop_sbatch_main($cmd, "footPeak", \@origFile, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir);
#sbatch_these($cmd, "footPeak", \@origFile, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir);

LOG($outLog, date() . "${LGN}SUCCESS$N: footPeak ran successfully!\n\n");


###############
# SUBROUTINES #
###############
sub record_options {
   my ($defOpts, $usrOpts, $usrOpts2, $other, $outLog, $logFile, $date, $uuid, $version_small) = @_;
   my $optPrint = "$0";
	my $optShort = "$0";
   foreach my $opt (sort keys %{$defOpts}) {
		next if $opt !~ /^.$/;
      my $val = $defOpts->{$opt};
		if (not defined $val) {print "opt=$opt, val undefnied!\n"}
      $optPrint .= $val eq "FALSE" ? "" : $val eq "" ? " -$opt " : $val eq "MYTRUE" ? " -$opt" : " -$opt $val";
		$optShort .= (not defined $usrOpts->{$opt}) ? "" : $val eq "FALSE" ? "" : $val eq "" ? "" : $val eq "MYTRUE" ? " -$opt" : " -$opt $val";
   }

   my  $param = "

$YW<-----------------------------------------------$N
${YW}Initializing...$N
>footPeak.pl version : $version_small
>Run Params
Date                : $date
Run ID              : $uuid
Run script short    : $optShort
Run script full     : $optPrint
>Options:
-d minDis  : $defOpts->{d}
-s grpSize : $defOpts->{s}
-k Dist    : $defOpts->{k}
-t thrshld : $defOpts->{t}
-w window  : $defOpts->{w}

>Run Params from footLoop.pl logfile=$logFile
footLoop Date        : $other->{date}
footLoop Run ID      : $other->{uuid}
footLoop Run script  : $other->{footLoop_runscript}
footLoop md5         : $other->{md5}
footLoop origDir     : $other->{origDir}
footLoop seqFile     : $defOpts->{seqFile}
";

if (defined $defOpts->{samFile}) {
	$param .= "
footLoop samFile     : $defOpts->{BAMFile}
footLoop samFixed    : $defOpts->{BAMFixed}
footLoop samFixedMD5 : $defOpts->{BAMFixedMD5}
";
}

if (defined $defOpts->{BAMFile}) {
	$param .= "
footLoop BAMFile     : $defOpts->{BAMFile}
footLoop BAMFixed    : $defOpts->{BAMFixed}
footLoop BAMFixedMD5 : $defOpts->{BAMFixedMD5}
";
}

	$param .= ">Options from footLoop.pl logfile=$logFile\n";
	my $optcount = 0;
	foreach my $opt (sort keys %{$defOpts}) {
		next if $opt !~ /^.$/;
		next if defined $usrOpts2->{$opt};
		#next if not defined $logs->{$opt}; next if $opt eq "origDir";
		$optcount ++;
		$param .= "-$opt : $defOpts->{$opt}\n" if defined $defOpts->{$opt};
	}
	$param .= "\n\n$YW----------------------------------------------->$N\n\n";

	LOG($outLog, $param);
}

sub parse_footLoop_logFile {
	my ($logFile, $date, $uuid, $footLoopFolder, $version_small) = @_;
	my @line = `cat $logFile`;
	my $paramsFile = "$footLoopFolder/.PARAMS";
	my ($defOpts, $usrOpts, $usrOpts2) = set_default_opts();
	my $inputFolder = $defOpts->{n}; 

	my @parline = `cat $paramsFile`;
	foreach my $parline (@parline) {
		if ($parline =~ /footLoop.pl,geneIndexFile,/) {
			($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,geneIndexFile,/) {
			($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
			($defOpts->{i}) = $parline =~ /geneIndexFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,seqFile,/) {
			($defOpts->{seqFile}) = $parline =~ /seqFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFile,/) {
			($defOpts->{samFile}) = $parline =~ /samFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,BAMFile,/) {
			($defOpts->{BAMFile}) = $parline =~ /BAMFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFixedFile,/) {
			($defOpts->{samFixed}) = $parline =~ /samFixedFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,BAMFixedFile,/) {
			($defOpts->{BAMFixed}) = $parline =~ /BAMFixedFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFixedFileMD5,/) {
			($defOpts->{samFixedMD5}) = $parline =~ /samFixedFileMD5,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,BAMFixedFileMD5,/) {
			($defOpts->{BAMFixedMD5}) = $parline =~ /BAMFixedFileMD5,(.+)$/;
		}
	}
	# %others contain other parameters from footLoop that aren't options (e.g. uuid)
	my ($other, $outLog); 
	
	foreach my $line (@line[0..@line-1]) {
		if ($line =~ /^\s*Date\s*:/) {
			($other->{date}) = $line =~ /^Date\s*:\s*([a-zA-Z0-9\:\-].+)$/;
			print "Undefined date from $line\n" and die if not defined $other->{date};
		}
		if ($line =~ /^\s*Run ID\s*:/) {
			($other->{uuid}) = $line =~ /^Run ID[ \t]+:[ \t]+([a-zA-Z0-9]+.+)$/;
			print "Undefined uuid from $line\n" and die if not defined $other->{uuid};
		}
		if ($line =~ /^\s*Run script\s*:/) {
			($other->{footLoop_runscript}) = $line =~ /^Run script[ ]+:(.+)$/;
			print "Undefined runscript from $line\n" and die if not defined $other->{footLoop_runscript};
			($defOpts, $other->{runscript}) = parse_runscript($defOpts, $usrOpts, $other->{footLoop_runscript});
		}
		if ($line =~ /^\s*Output\s*:/) {
			my ($value) = $line =~ /Output.+(\.0_orig_\w{32})/;
			die "Died at line=$line, value=?\n" if not defined $value;
			$value = $footLoopFolder . "/$value";
			die "Died at line=$line, value=?\n" if not -d $value;
			($other->{origDir}) = $value;
			print "Undefined origDir from input=$inputFolder, line=$line\n" and die if not defined $other->{origDir};
			print "origDir $other->{origDir} does not exist!\n" and die if not -d $other->{origDir};
			$defOpts->{origDir} = $other->{origDir};
			($other->{md5}) = $other->{origDir} =~ /\.0_orig_(\w{32})/;
			print "Can't parse md5 from outdir (outdir=$defOpts->{origDir})\n" and die if not defined $other->{md5};
		}
#		if ($line =~ /geneIndexFile=/) {
#			($defOpts->{geneIndexFile}) = $line =~ /geneIndexFile=(.+)$/ if $line !~ /,gene=.+,beg=\d+,end=\d+$/;
#			($defOpts->{geneIndexFile}) = $line =~ /geneIndexFile=(.+),gene=.+,beg=\d+,end=\d+$/ if $line =~ /,gene=.+,beg=\d+,end=\d+$/;
#			$defOpts->{geneIndexFile} = $footLoopFolder . "/" .  getFilename($defOpts->{geneIndexFile});
#		}
		if ($line =~ /^!\w+=/) {
			my ($param, $value) = $line =~ /^!(\w+)=(.+)$/;
			my $param2 = defined $param ? $param : "__UNDEF__";
			my $value2 = defined $value ? $value : "__UNDEF__";
			if ($value =~ /\//) {
				if ($value =~ /\/?\.0_orig\w{32}/) {
					($value) = $value =~ /^.+(\/?\.0_orig_\w{32})/; 
					$value = $footLoopFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
				if ($value =~ /\/\.geneIndex/) {
					($value) = $value =~ /^.+(\/\.geneIndex.+)$/; 
					$value = $footLoopFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
				else {
					($value) = getFilename($value, 'full');
					$value = $footLoopFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
			}
			print "$param = $value\n";
			print "Cannot parse param=$param2 and value=$value2 from line=$line\n" and die if not defined $param or not defined $value;
			print "$param file $value does not exist!\n" and die if $value =~ /\/+/ and not -e $value;
			($defOpts->{$param}) = $value if $param ne "n";
		}
	}
	if (not defined $other->{origDir}) {
		($other->{origDir}) = $defOpts->{n};
		print "Undefined origDir from input=$inputFolder\n" and die if not defined $other->{origDir};
		print "origDir $other->{origDir} does not exist!\n" and die if not -d $other->{origDir};
		$defOpts->{origDir} = $other->{origDir};
		($other->{md5}) = "NA"; #$other->{origDir} =~ /\.0_orig_(\w{32})/;
	}
	$defOpts->{o} = $defOpts->{n} if not defined $opt_o;
	makedir($defOpts->{o}) if not -d $defOpts->{o};
	#die "opt = $opt_o = $defOpts->{o}\n";
	open ($outLog, ">", "$defOpts->{o}/footPeak_logFile.txt") or print "Failed to write to $defOpts->{o}/footPeak_logFile.txt: $!\n" and die;
	if (defined $opt_O) {
		$defOpts->{origDir} = $opt_O;
	}
	elsif (not defined $other->{origDir}) {
		$defOpts->{origDir} = $defOpts->{o};
	}
	elsif (not defined $defOpts->{origDir}) {
		$defOpts->{origDir} = $defOpts->{o};
	}
	$other->{md5} = `which md5sum`;
	chomp($other->{md5});

	if (defined $defOpts->{BAMFile}) {
		$defOpts->{BAMFixed} = $defOpts->{BAMFile} . ".fixed.gz";
		$defOpts->{BAMFixedMD5} = getMD5($defOpts->{BAMFile} . ".fixed.gz");
	}
	
	record_options($defOpts, $usrOpts, $usrOpts2, $other, $outLog, $logFile, $date, $uuid, $version_small);

#	print "Output = $defOpts->{o}\n";
	return($defOpts, $outLog);
}

sub parse_runscript {
	my ($defOpts, $usrOpts, $runscripts) = @_;
	my @runscripts = split(" ", $runscripts);
	my $log;
	for (my $i = 0; $i < @runscripts; $i++) {
		my ($opt, $val);
		if ($runscripts[$i] =~ /^\-[a-zA-Z]$/) {
			($opt) = $runscripts[$i] =~ /^\-(.)$/;
			if ($i < @runscripts-1) {
				$val = $runscripts[$i+1];
			}
			else {
				$val = "MYTRUE";
			}
			if (not defined $defOpts->{$opt}) {
				print "parse_runscript: opt $opt in footLoop logFile.txt Run Script does not exist in footPeak getopt options!\n" and next;
			}
		}
		elsif ($i == 0) {
			next;#print "script = $runscripts[0]\n" and next;
		}
		else {
			print "parse_runscript: opt \$opt in footLoop logFile.txt Run Script does not look like a getopt options (not a '-a')!\n" and next;
		}
		if (defined $val and $val =~ /^\-[a-zA-Z]$/) {
			my ($val1) = $val =~ /^\-(.)$/;
			if (defined $defOpts->{$val1}) {
				$val = "MYTRUE";
			}
			else {
				$i ++;
			}
		}
		else {
			$i ++;
		}
		my $opt2 = defined $opt ? $opt : "__UNDEF__";
		my $val2 = defined $val ? $val : "__UNDEF__";
		print "i=$i Undefined opt=$opt2 val=$val2 line=$runscripts[$i]\n" and die if not defined $val or not defined $opt;
		next if $opt eq "n";
		if ($opt eq "l") {
			$defOpts->{label} = $val;
			$log->{label} = 1;
			next;
		}
		$defOpts->{$opt} = $val;# eq "MYTRUE" ? "MYTRUE" : $val;
		$log->{$opt} = 1;
	}
	my $runscript = "$0";
	foreach my $opts (sort keys %{$defOpts}) {
		next if $opts =~ /^[A-Z]$/;
		next if not defined $usrOpts->{$opts};
		next if $defOpts->{$opts} eq "FALSE";
		$runscript .= " -$opts $defOpts->{$opts}";# if defined $usrOpts->{$opts};
	}
	foreach my $opts (sort keys %{$defOpts}) {
		next if $opts !~ /^[A-Z]$/;
		next if not defined $usrOpts->{$opts};
		next if $defOpts->{$opts} eq "FALSE";
		$runscript .= " -$opts $defOpts->{$opts}";# if defined $usrOpts->{$opts};
	}
	foreach my $opts (sort keys %{$defOpts}) {
		next if $opts =~ /^[A-Z]$/;
		next if defined $usrOpts->{$opts};
		next if $defOpts->{$opts} eq "FALSE";
		$runscript .= " -$opts $defOpts->{$opts}";# if defined $usrOpts->{$opts};
	}
	foreach my $opts (sort keys %{$defOpts}) {
		next if $opts !~ /^[A-Z]$/;
		next if defined $usrOpts->{$opts};
		next if $defOpts->{$opts} eq "FALSE";
		$runscript .= " -$opts $defOpts->{$opts}";# if defined $usrOpts->{$opts};
	}
#	print $runscript and die;
#	print "$defOpts->{i}\n";die;
	return ($defOpts, $runscript);
}

sub set_default_opts {
	# To define default values as well as record purposes, we create 4 hashes, %defOpts, %usrOpts, %other, and %log.
	# %defOpts contains options that came from footLoop logfile.txt and came from user inputs from this script (footPeak.pl)
	# %defOpts stores default values, %usrOpts stores user inputs.
	my %defOpts =
		(
		'd' => '250'     ,'g' => ''        ,'i' => ''        ,'k' => '50'      ,
		'l' => '100'     ,'n' => ''        ,'q' => '0'       ,'r' => ''        ,
		's' => '200'     ,'t' => '55'      ,'w' => '20'      ,'x' => '0'       ,
		'y' => '0'       ,'K' => '2'       ,'L' => ''       ,'A' => ''         ,
		'o' => 'RESULTS' ,'G' => ''        ,'Z' => ''       ,'J' => ''         ,
		'p' => ''
	);


	# %usrOpts contains options only from this script (footPeak.pl) 
	# For each undefined values in %usrOpts (user didn't input anything), we use default value in %defOpts.
	my %usrOpts =
		(
		'd' => $opt_d    ,'g' => $opt_g    ,'h' => 'FALSE'   ,'k' => $opt_k    ,
		'l' => $opt_l    ,'n' => $opt_n    ,'p' => $opt_p    ,'s' => $opt_s    ,
		't' => $opt_t    ,'w' => $opt_w    ,'K' => $opt_K    ,'A' => 'FALSE'   ,
		'o' => $opt_o    ,'G' => 'FALSE'   ,'L' => ''
		);

	my %usrOpts2 =
		(
		'd' => 'd'    ,'g' => 'g'    ,'h' => 'FALSE'    ,'k' => 'k'    ,
		'l' => 'l'    ,'n' => 'n'    ,'p' => 'p'    		,'s' => 's'    ,
		't' => 't'    ,'w' => 'w'    ,'K' => 'K'        ,'G' => 'G'    ,
      'A' => 'A'    ,'L' => 'L'
		);

	print_default_opts(\%defOpts, \%usrOpts) and die if @ARGV == 1 and $ARGV[0] eq "ex";

	# set opt in footPeak's %defOpts to be what user inputted based on %usrOpts.
	# if doesn't exist, then use default value in %defOpts.
	foreach my $opt (sort keys %usrOpts) {
		next if not defined $usrOpts{$opt};
		$defOpts{$opt} = $usrOpts{$opt};
		if (-d $usrOpts{$opt}) {
			$usrOpts{$opt} = "./$usrOpts{$opt}" if $usrOpts{$opt} !~ /\/.+$/;
			$defOpts{$opt} = getFullpath($usrOpts{$opt});#, 1) . "/" if $opt eq "o";
		}
	}
	# This below is to print
	return(\%defOpts, \%usrOpts, \%usrOpts2);
}

sub print_default_opts {
	my ($defOpts, $usrOpts) = @_;
	my ($defOptsCount, $usrOptsCount) = (0,0);
	my $defOptsPrint = "\tmy \%defOpts =\n\t\t(";
	foreach my $opt (sort keys %{$defOpts}) {
		next if $opt =~ /\-?[A-Z]$/;
		my $val = $defOpts->{$opt} =~ /^\-\d+$/ ? $defOpts->{$opt} : "\'$defOpts->{$opt}\'";
		   $val = $val . join("", (" ") x (10 - length($val))) . ",";
		$defOptsPrint .= "\n\t\t" if $defOptsCount % 4 == 0;
		$defOptsPrint .= "\'$opt\' => $val";
		$defOptsCount ++;
	}
	foreach my $opt (sort keys %{$defOpts}) {
		next if $opt !~ /\-?[A-Z]$/;
		my $val = $defOpts->{$opt} =~ /^\-\d+$/ ? $defOpts->{$opt} : "\'$defOpts->{$opt}\'";
		   $val = $val . join("", (" ") x (10 - length($val))) . ",";
		$defOptsPrint .= "\n\t\t" if $defOptsCount % 4 == 0;
		$defOptsPrint .= "\'$opt\' => $val";
		$defOptsCount ++;
	}
	$defOptsPrint =~ s/,$/\n\t);/;

	print "$defOptsPrint\n\n";
	parse_getopt();
}

sub parse_getopt {
	my ($usrOpts) = `grep 'getopt' $0`; 
	die "Undef usropts1\n" if not defined $usrOpts;
	chomp($usrOpts); 
	my ($usrOpts2) = $usrOpts =~ /^.+opts\(\"(.+)\"\);/; 
	die "Undef usropts2 (usrOpts = $usrOpts)\n" if not defined $usrOpts2;
	my @usrOpts = split(":", $usrOpts2);
	my %usrOpts;
	my $print = "\tmy \%usrOpts =\n\t\t(";
	for (my $i = 0; $i < @usrOpts; $i++) {
		$usrOpts[$i] =~ s/:$//;
		my @opts = split("", $usrOpts[$i]);
		my $val = "\'FALSE\'";
		for (my $j = 0; $j < @opts-1; $j++) {
			my $opt = $opts[$j];
			next if $opt eq "v";
			$usrOpts{$opt} .= "'$opt' => $val" . join("", (" ") x (10 - length($val))) . ",";
		}
		$val = "\$opt_$opts[@opts-1]";
		$usrOpts{$opts[@opts-1]} .= "'$opts[@opts-1]' => $val" . join("", (" ") x (10 - length($val))) . ",";
	}
	my $count = 0;
	foreach my $opt (sort keys %usrOpts) {
		next if $opt =~ /^[A-Z]$/;
		$print .= "\n\t\t" if $count % 4 == 0;
		$print .= $usrOpts{$opt};
		$count ++;
	}
	foreach my $opt (sort keys %usrOpts) {
		next if $opt !~ /^[A-Z]$/;
		$print .= "\n\t\t" if $count % 4 == 0;
		$print .= $usrOpts{$opt};
		$count ++;
	}
	$print =~ s/,$/\n\t\t);\n/;
	
	print "$print";
	die;
}
sub parse_indexFile_and_seqFile {
	my ($indexFile, $seqFile, $outLog) = @_;
	my $SEQ;

	my @indexLine = `cat $indexFile`;
	my $totalindexLine = @indexLine;
	LOG($outLog, "Index File = $indexFile\n");
	foreach my $line (@indexLine) {
		chomp($line);
		my ($chr, $beg, $end, $def, $zero, $strand) = split("\t", $line);
		$def = uc($def);
		$SEQ->{$def}{coor} = "$chr\t$beg\t$end\t$def\t$zero\t$strand";
#		LOG($outLog, "indexFile:\tSEQ -> {gene=$def}{coor} is defined!\n");
		LOG($outLog, "def=$def, coor=$chr, $beg, $end, $def, $zero, $strand\n","NA");
	}
	open (my $in, "<", $seqFile) or DIELOG($outLog, "Failed to read from $LCY$seqFile$N: $!\n");
	my $linecount = 0;
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def; $def =~ s/^>//; $def = uc($def);
		my @seq = split("", $entry->seq);
		$SEQ->{$def}{seq} = \@seq;
		$SEQ->{$def}{loc} = findCGPos(\@seq);
		LOG($outLog, "seqFile:\tSEQ -> {gene=$def}{seq} and {loc} is defined!\n","NA");
		$linecount ++;
	}
	close $in;
	LOG($outLog, "\n" . date() . "${LPR}parse_indexFile_and_seqFile$N: Parsed $LGN$totalindexLine$N genes from .bed and $LGN$linecount$N references from .fa!\n");
	return $SEQ;
}

sub findCGPos {
   my $seq = $_[0];
   my $data;
   for (my $i = 0; $i < @{$seq}; $i++) {
      if ($seq->[$i] eq "C") {
         my $pos = defined $data->{pos1} ? (keys %{$data->{pos1}}) : 0;
         $data->{pos1}{$pos} = $i;
         $data->{pos2}{$i} = $pos;
      }
      if ($seq->[$i] eq "G") {
         my $neg = defined $data->{neg1} ? (keys %{$data->{neg1}}) : 0;
         $data->{neg1}{$neg} = $i;
         $data->{neg2}{$i} = $neg;
      }
   }
   return $data;
}


sub check_sanity {
	my ($footLoopFolder) = @_;
	my ($usage, $usage_long) = usage();
	if (defined $opt_h or not defined $footLoopFolder or not -d $footLoopFolder) {
		print $usage;
		print $usage_long if defined $opt_h;
		print "${LRD}ERROR:$N Please define -n (output dir from footLoop.pl)\n" if not defined $footLoopFolder or not -d $footLoopFolder;
		print "${LRD}ERROR:$N Please define -o (output dir of footPeak.pl)\n" if not defined $opt_o;
		die "${YW}----------------------$N\n\n";
	}
	LOG($outLog, "Please run footloop.pl first! ($footLoopFolder/logFile.txt does not exists)\n") and die if not -e "$footLoopFolder/logFile.txt";

	return $usage;
}

# <-------------------------------

sub usage {

my ($scriptFolder, $scriptName) = getFilename($0, "folderfull");
####### Usage #######
	my $usage = "
${YW}----------------------$N
$YW|$N $scriptName $version_small  $YW|$N
${YW}----------------------$N

Usage: $YW$scriptName$LGN [Options: -w|-t|-G]$N -n $LCY<footLoop output folder>$N -o $LCY<output png dir>$N

Options$LGN [default]$N:
-w: window size$LGN [20]$N
-t: threshold in \%$LGN [55]$N
-G: only process this gene$LGN [NA]$N
-l: min length of peak$LGN [100]$N

footLoop folder: $scriptFolder$N

${YW}---------------------$N
";

####### Usage Long #######
my $usage_long = "
${YW}Long Usage$N:

-d of 100:
Say there's another peak >> WITHIN SAME READ << that begins at position 280 and ending at position 500.
This peak is 'too close' to the first peak as its beginning (280) is less than 100bp away from the first peak's end (255)
Therefore these two peak -will be merged-
New peak begin at 108 end at 500

-s of 200:
	beg of this peak = integer of (108/200) = 0
	end of this peak = integer of (255/200) = 1

The distance of this start/end bin with another peak's start/end bin determines its cluster.
This means this peak will be -very- close with other peaks that also begin at bin 0 and end at bin 1
It will be slightly further with peak begin at bin 0 and end at bin 2
It will be weakly grouped with peak beginning at bin 0 and end at bin 99

-k of 50:
beg = 108-50 to 108+50 = 58 to 158
mid = 108+50 to 255-50 = 158 to 205
end = 255-50 to 255+50 = 205 to 305

${YW}---------------------$N

";

##########################

	return($usage, $usage_long);
}

# -------------------------------->

# 0 is bad or not data (was 6)
# 1 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak



#my ($i, @origFile, $out, $outDir, $window, $threshold, $genewant, $seqFile, $totalorigFile, $label);
#sbatch_these($bismark_cmd, "bismark", "bam", \@fastqFiles, $max_slurm_job, $outLog);
=comment
for (my $i = 0; $i < @origFile; $i++) {
	my $peakFile = $origFile[$i];
	LOG($outLog, date() . "$LGN$i$N/$LGN$totalorigFile$N:$LPR footPeak.pl processing$N $YW$peakFile$N\n");

	if (defined $opt_G and $peakFile !~ /$opt_G/i) {
		LOG($outLog, date() . "-> Skipped $LCY$peakFile$N as it doesn't contain $LGN-G $opt_G$N\n","NA");
		next;
	}

	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
	$peakFilename =~ s/.filtered.gz$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	DIELOG($outLog, "gene=$gene seq->gene{seq} isn't defined!\n") if not defined $SEQ->{$gene} or not defined $SEQ->{$gene}{seq};
	my $cmd = "footPeak_sbatch_2.pl $indexFile $seqFile $i FILENAME $outDir $window $threshold $totalorigFile $label $genewant $minDis $minLen $version_small";
	sbatch_these($cmd, "");
	LOG($outLog, $cmdout);
	
	print "DEBUG! last\n" if defined $opt_0;
	last if defined $opt_0;
}
LOG($outLog, date() . "${LGN}SUCCESS$N: footPeak ran successfully!\n\n");
=cut
