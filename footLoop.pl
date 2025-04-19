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

use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars   qw($opt_v $opt_r $opt_g $opt_G $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f $opt_l $opt_e $opt_z $opt_9 $opt_o $opt_J $opt_0 $opt_a $opt_8 $opt_f);
my @opts = qw(       $opt_r $opt_g $opt_i $opt_G $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_l $opt_e $opt_z $opt_9 $opt_o $opt_J $opt_0 $opt_a $opt_8 $opt_f);
getopts("vr:g:G:i:n:L:x:y:q:HhZF:p:l:ez9o:J:0a8f");

BEGIN {
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";
}
use myFootLib;
use FAlite;

##################
# 0. Check Sanity #
##################

my %OPTS  = ('r' => $opt_r, 'g' => $opt_g, 'i' => $opt_i, 'n' => $opt_n, 'G' => $opt_G,
				 'L' => $opt_L, 'x' => $opt_x, 'y' => $opt_y, 'p' => $opt_p, 'J' => $opt_J,
				 'q' => $opt_q, 'F' => 'NONE', 'Z' => 'NONE', 'o' => $opt_o, '0' => $opt_0,
				 'p' => 'NONE', 'q' => $opt_q, 'l' => $opt_l);
my %OPTS2 = ('p' => $opt_p, 'Z' => $opt_Z, 'F' => $opt_F, 'a' => $opt_a);

my $footLoop_script_folder = dirname(dirname abs_path $0) . "/footLoop/";
my $fixBAMFile_script = "$footLoop_script_folder/bin/footLoop_2_fixBAMFile.pl";
my $filterBAMFile_script = "$footLoop_script_folder/bin/footLoop_3_filterBAMFile.pl";

my $version = "unk";
my $md5script = "unk";
check_sanity(\%OPTS, \%OPTS2);
($footLoop_script_folder, $version, $md5script) = check_software();

print "version $YW$version$N\n";
$opt_F = "1,2,3,4,5,6,7,8,9" if defined $opt_f;
###################
# 1. Define Input #
###################
my $snakemake;
my $ambig = $opt_a;
my ($readFile, $genomeFile, $geneIndexFile, $outDir) = getFullpathAll($opt_r, $opt_g, $opt_i, $opt_n);
my ($readFilename)  = getFilename($opt_r, "full");
my ($geneIndexName) = getFilename($geneIndexFile);
my $minMapQ    = (not defined($opt_q)) ? 0  : $opt_q;
$opt_L = "95p" if not defined $opt_L;
my ($minReadL)   = $opt_L =~ /p$/i ? $opt_L =~ /^(.+)p$/i : $opt_L;
my $bufferL    = (not defined($opt_x)) ? 0  : $opt_x;
my $bufferR    = (not defined($opt_y)) ? 0  : $opt_y;
my $bufferLen  = $bufferR - $bufferL;
my $readName   = getFilename($readFile, "full");
my $bowtieOpt = (defined $opt_Z) ? "--bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8" : "--bowtie2 --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3";
   $bowtieOpt = (defined $opt_z) ? "--bowtie2 --rdg 3,2 --rfg 3,2 --score_min L,0,-0.5" : $bowtieOpt;
my $BAMFile	   = ($readFile =~ /.f(ast)?q(.gz)?$/i) ? $outDir .  "/$readName\_bismark_bt2.bam" : $readFile;
my $uuid       = getuuid();
my $date       = getDate();
my $filteredDir    = $outDir . "/.0_filteredBAM";
my $parallel    = defined $opt_p ? $opt_p : 1000;
my $max_slurm_job = defined $opt_J ? $opt_J : 500;
$max_slurm_job = 500 if $max_slurm_job < 0;

my %force;
if (defined $opt_F) {
	my @forces = split(",", $opt_F);
	foreach my $force (@forces) {
		$force{$force} = 1;
	}
	if (defined $force{1}) {print "1: rerun all fasta preprocessin\n";}
	if (defined $force{2}) {print "2: delete all part.gz\n";}
	if (defined $force{3}) {print "3: split Fastq into part.gz\n";}
	if (defined $force{4}) {print "4: bismark part.gz\n";}
	if (defined $force{5}) {print "5: filter bam ($fixBAMFile_script)\n";}
	if (defined $force{6}) {print "6: log bam ($filterBAMFile_script)\n";}
	if (defined $force{7}) {print "7: merge fixed.gz\n";}
	if (defined $force{8}) {print "8: merge filtered.gz\n";}
	if (defined $force{9}) {print "9: merge bam\n";}
}

# Make directory
makedir($outDir);

# Make log file
my $logFile = "$outDir/logFile.txt";
my $bigCMDFile = "$outDir/bigCMD.sh";
open(my $outReadLog, ">", "$outDir/.PARAMS") or print "\n${LRD}FATAL!$N Failed to write to $LCY$outDir/.PARAMS$N: $!\n" and die;
open(my $outLog, '>', $logFile) or print "\n${LRD}FATAL!$N Failed to write to $LCY$logFile$N: $!\n" and die;
open(my $outBigCMD, '>', $bigCMDFile) or print "\n${LRD}FATAL!$N Failed to write to $LCY$bigCMDFile$N: $!\n" and die;
LOG($outReadLog, "footLoop.pl,uuid,$uuid\n");
LOG($outReadLog, "footLoop.pl,date,$date\n");
# Record all options

record_options(\%OPTS, \%OPTS2, $outReadLog, $version, $outLog);
#######################################
# 1. Preprocess Index and Fasta Files #
#######################################

# Add buffers to geneIndexFile ($footLoop_script_folder/lib/bedtools_bed_change.pl) and get their sequences using fastaFromBed, then get their info from fasta=$seqFile

my ($SEQ, $geneIndexHash, $seqFile, $bismark_geneIndexDir, $geneIndexFile2) = parse_geneIndexFile($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog);

#####################################
# 2. Run Bismark Genome Preparation #
#####################################

# Make Bismark Index
my $forcerun = defined $force{1} ? 1 : 0;
bismark_genome_preparation($seqFile, $bismark_geneIndexDir, $bowtieOpt, $forcerun, $outLog);

##################
# 3. Run Bismark #
##################

# Run Bismark
if (not defined $opt_p) {
	($BAMFile) = run_bismark($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile2, $SEQ, $outReadLog, $bismark_geneIndexDir, $bowtieOpt, $outLog);
}
else {
	($BAMFile) = run_bismark_parallel($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile2, $SEQ, $outReadLog, $bismark_geneIndexDir, $bowtieOpt, $outLog);
}
my ($BAMFileName) = getFilename($BAMFile);

LOG($outLog, "\n" . date() . "${LGN}SUCCESS!!$N: footLoop.pl ran successfuly\n");

###############
# SUBROUTINES #
###############


sub bismark_genome_preparation {
	my ($geneIndexFa, $bismark_geneIndexDir, $bowtieOpt, $forcerun, $outLog) = @_;

	my $bismark_genome_preparation_cmd = "bismark_genome_preparation --bowtie2 $bismark_geneIndexDir > $bismark_geneIndexDir/LOG.txt 2>&1";

	LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: Running bismark_genome_preparation\n");
	my ($geneIndexFaMD5, $temp, $geneIndexFaTempMD5File)  = getMD5($geneIndexFa);
	my $oldMD5File = "$bismark_geneIndexDir/Bisulfite_Genome/MD5SUM";

	my $run_boolean = "\t${LGN}WAS NOT RUN (Using previously made files)$N:";
	my ($check, $md5sum, $md5sum2) = (0,$geneIndexFaMD5,"md5sum2");
	my $bismark_geneIndexDir_exist = (-d "$bismark_geneIndexDir/Bisulfite_Genome/" and -e $geneIndexFaTempMD5File) ? 1 : 0;


	if ($bismark_geneIndexDir_exist == 1 and $forcerun eq 0) { # in case in the future we implement -G for bismark folder
		LOG($outLog, "\t#Older bismark folder Bisulfite_Genome $CY$bismark_geneIndexDir/Bisulfite_Genome/$N exist! Checking MD5 if they're the same as current fasta file.\n","NA");
		($md5sum)  = `cat $oldMD5File` =~ /^(\w+)($|[\t ].+)$/;
		($md5sum2) = getMD5($geneIndexFa);
		$bismark_geneIndexDir_exist = 0 if not defined $md5sum or not defined $md5sum2;
		if (defined $md5sum and defined $md5sum2) {
			$bismark_geneIndexDir_exist = 0 if $md5sum ne $md5sum2;
		}
	}
	$run_boolean = "\t${LGN}RAN$N:" if $bismark_geneIndexDir_exist == 0 or $forcerun eq 1;

	if ($bismark_geneIndexDir_exist == 0 or $forcerun eq 1) {
		print_cmd($bismark_genome_preparation_cmd, $outLog);
		LOG($outLog, "\t#Either bismark folder didn't exist or older bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n") if defined $md5sum and $bismark_geneIndexDir_exist == 0;
		LOG($outLog, "\t#Rerunning bismark_genome_preparation due to user request (-F 1)\n") if $forcerun eq 1;
		system($bismark_genome_preparation_cmd) == 0 or die "Failed to run bismark genome preparation: $!\n";
		LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: ${GN}SUCCESS$N: $CY$bismark_geneIndexDir\/Bisulfite_Genome$N made\n");
	}
	else { 
		LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: ${GN}SUCCESS$N: $CY$bismark_geneIndexDir\/Bisulfite_Genome$N already exist (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n");
	}
	system("echo $geneIndexFaMD5 > $bismark_geneIndexDir/Bisulfite_Genome/MD5SUM") == 0 or die "Failed to generate $LCY$bismark_geneIndexDir/Bisulfite_Genome/MD5SUM/$N: $!\n";
}

sub run_bismark {
	my ($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile, $SEQ, $outReadLog, $bismark_geneIndexDir, $bowtieOpt, $outLog) = @_;
	my $mapLog = "";
	$BAMFile =~ s/(.fq|.fastq|.fq.gz|.fastq.gz)_bismark_bt2/_bismark_bt2/;
	my $ext = "bam";
	my ($BAMFilename) = getFilename($BAMFile, "full");
	
	my $run_boolean = "\n\t${LGN}WAS NOT RUN$N:${YW} ";

	$run_boolean = "\n$YW\t ";

	my $outFolder = $outDir;
	LOG($outLog, "\t#Not parallel\n");

	LOG($outLog, "==================================================\n");			

	###########
	# BISMARK #
	###########

	my $bismark_cmd       = "bismark --non_directional -o $outFolder --temp_dir $outFolder --ambig_bam --ambiguous --unmapped $bowtieOpt $bismark_geneIndexDir $readFile > $outFolder/.bismark_log 2>&1";
	my $bismark_cmd_print = "bismark --non_directional -o \$outFolder --temp_dir \$outFolder --ambig_bam --ambiguous --unmapped \$bowtieOpt \$bismark_geneIndexDir \$readFile > \$outFolder/.bismark_log 2>&1";
	$run_boolean .= " ::: $bismark_cmd :::";
	my $totalbamFiles = 1;
	my $BAMFileDone = "$BAMFile.Done";
	if (not -e $BAMFileDone or defined $force{4}) {
		LOG($outLog, "$LCY\n\t$bismark_cmd\n$N\n");
		LOG($outReadLog, "footLoop.pl,bismark,$bismark_cmd\n","NA");
		print_cmd($bismark_cmd, $outLog);
		my $result = system($bismark_cmd);

		if ($result != 0) {
			LOG($outLog, "\t\t${LRD}bismark failed!$N In case bisulfte_genome is corrupted, we're re-running bismark_genome_preparation:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_geneIndexDir\n");
			LOG($outReadLog, "$LCY\n\tfootLoop.pl,bismark,bismark_genome_preparation --bowtie2 $bismark_geneIndexDir\n$N","NA");
			bismark_genome_preparation($seqFile, $bismark_geneIndexDir, $bowtieOpt, $forcerun, $outLog);
			LOG($outReadLog, "footLoop.pl,bismark,$bismark_cmd");
			system($bismark_cmd) == 0 or DIELOG($outLog, "$LRD!!!$N\tFailed to run bismark: $!\n");
		}
		LOG($outLog, "\n" . date() . "${LPR}bismark .fastq.gz .bam$N: ${GN}SUCCESS$N: Created $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$BAMFile\n");
		system("touch $BAMFileDone") if not -e $BAMFileDone;
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}bismark .fastq.gz .bam$N: ${GN}SUCCESS$N: Using previously run $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$BAMFile$N\n");
	}

	LOG($outLog, "==================================================\n");			
	################
	# FIX BAM FILE #
	################
	LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: Running ${LCY}$fixBAMFile_script$N\n");
	my $fixedBAMFile = "$BAMFile.fixed.gz";
	my $fixedBAMFileDone = "$BAMFile\_fixBAM.done";
	my $fix_BAMFile_cmd = "$fixBAMFile_script -b $BAMFile -g $seqFile -i $geneIndexFile -o $outFolder";
	my $totalfixedFiles = 1;
	if (not -e $fixedBAMFileDone or defined $force{5}) {
		print_cmd("$fix_BAMFile_cmd", $outLog);
		system($fix_BAMFile_cmd) == 0 or DIELOG($outLog, "footLoop.pl::run_bismark: Failed to run fixBAMFile: $LRD$!$N\ncmd=$LCY$fix_BAMFile_cmd$N\n\n");
		LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: ${GN}SUCCESS$N: Created $LGN$totalfixedFiles$N .fixed.gz fixed Bam Files in $LCY$fixedBAMFile$N\n");
		system("touch $fixedBAMFileDone") if not -e $fixedBAMFileDone;
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: ${GN}SUCCESS$N: Using previously run $LGN$totalfixedFiles$N .fixed.gz fixed Bam Files in $LCY$fixedBAMFile$N\n");
	}
	LOG($outLog, "==================================================\n");			
	
	###################
	# FILTER BAM FILE #
	###################
	LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: Running ${LCY}$filterBAMFile_script$N\n");

	my $filteredBAMFile = "$BAMFile.log.txt";
	my $filteredBAMFileDone = "$BAMFile\_filterBAM.done";
	my $totalfilterBAMFile = 1;
	my $filter_BAMFile_cmd = "$filterBAMFile_script -r $opt_r -b $BAMFile -g $seqFile -L $opt_L -O $outFolder -o $outFolder -i $geneIndexFile2";
	   $filter_BAMFile_cmd .= " -q $opt_q" if defined $opt_q;
	   $filter_BAMFile_cmd .= " -x $opt_x" if defined $opt_x;
	   $filter_BAMFile_cmd .= " -y $opt_y" if defined $opt_y;


	if (defined $force{6} or not -e $filteredBAMFileDone) {
		print_cmd("$filter_BAMFile_cmd", $outLog);
		system("$filter_BAMFile_cmd") == 0 or DIELOG($outLog, "footLoop.pl::run_bismark: Failed to run filterBAMFile: $LRD$!$N\ncmd=$LCY$filter_BAMFile_cmd$N\n\n");
		LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: ${GN}SUCCESS$N: Created $LGN$totalfilterBAMFile$N *.log.txt files in $LCY$filteredBAMFile$N\n");
		system("touch $filteredBAMFileDone") if not -e $filteredBAMFileDone;
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: ${GN}SUCCESS$N: Using previously run $LGN$totalfilterBAMFile$N *.log.txt in $LCY$filteredBAMFile$N\n");
	}
	LOG($outLog, "==================================================\n");			
	##################

	#########################
	# MERGE BISMARK REPORTS #
	#########################
	LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: Merging reports in outFolder=$LCY$outFolder/*_SE_report.txt$N\n");
	print_cmd("#footLoop.pl::parse_bismark_report(\$part_bismark_report, \$outLog)", $outLog);
	my $report_not_found = 0;
	my $bamFiles = $BAMFile;
	my ($partFilename) = $bamFiles =~ /^(.+).bam$/;
	my $part_bismark_report = $partFilename . "_SE_report.txt";
	if (not -e $part_bismark_report) {
		$report_not_found ++;
		LOG($outLog, "\nfootLoop.pl::run_bismark: Can't find report file $LCY$part_bismark_report$N\n");
		next;
	}
	my ($header, $report) = parse_bismark_report($part_bismark_report, $outLog);
	my @header = @{$header};
	my @report = @{$report};
	my @HEADER = @header;
	my @REPORT;
	LOG($outLog, "\t$LCY$part_bismark_report$N:\n\t","NA");
	for (my $q = 0; $q < @header; $q++) {
		# just next coz we're logging, not that important
		next if not defined $header[$q];
		next if $header[$q] ne $HEADER[$q];
		next if not defined $report[0];
		next if not defined $report[$q];
		$REPORT[$q] += $report[$q] if $header[$q] !~ /^perc_/;
		$REPORT[$q] += $report[0]*$report[$q]/100 if $header[$q] =~ /^perc_/; 
		LOG($outLog, "q=$q, $report[0]*$report[$q]/100 = $REPORT[$q]\n","NA") if $header[$q] =~ /^perc_/;
	}
	LOG($outLog, "\n","NA");
	my $merged_report_to_log = "";
	my ($bismark_report) = $BAMFile =~ /^(.+).$ext/; $bismark_report .= "_footLoop_report.txt";
	my $bismark_report_exist = -e $bismark_report ? 1 : 0;

	open (my $outbismark_report, ">", $bismark_report) or DIELOG($outLog, "footLoop.pl::run_bismark: parallel: Failed to write to bismark report $LCY$bismark_report$N: $!\n");
	for (my $q = 0; $q < @HEADER; $q++) {
		$REPORT[$q] = int(10000*$REPORT[$q]/($REPORT[0])+0.5)/100 if $HEADER[$q] =~ /^perc_/; 
		print $outbismark_report "$HEADER[$q]\t$REPORT[$q]\n";
		$merged_report_to_log .= "$LCY$HEADER[$q]$N=$LGN$REPORT[$q]$N\n";
	}
	close $outbismark_report;

	$mapLog  = "footLoop.pl,map," . "header\t" . join("\t", @HEADER) . "\tfootLoop_outDir\tuuid\n";
   $mapLog .= "footLoop.pl,map," . "record\t" . join("\t", @REPORT) . "\t$outDir\t$uuid\n";
	if ($bismark_report_exist == 0) {
		LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: ${GN}SUCCESS$N: Merged $LGN$totalbamFiles$N bismark reports in $LCY$outFolder*/*_SE_report.txt$N, stats:\n$merged_report_to_log\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: ${GN}SUCCESS$N: Using previously merged $LGN$totalbamFiles$N bismark reports in $LCY$outFolder*/*_SE_report.txt$N, stats:\n\n$merged_report_to_log\n");
	}
	LOG($outLog, "\t#MAP_RESULT\n$LGN$mapLog$N\n");
	LOG($outLog, "==================================================\n");			

	# print to $outReadLog
	LOG($outReadLog, "footLoop.pl,BAMFile,$outDir/$BAMFilename\n","NA");
	LOG($outReadLog, $mapLog,"NA");

	my ($new1, $new2, $new3) = (newer($BAMFile,$fixedBAMFile),newer($BAMFile,$filteredBAMFile),newer($fixedBAMFile,$fixedBAMFile));
	my @Ftorun = ();
	if (newer($BAMFile,$fixedBAMFile) == 1) {
		LOG($outLog, "\n" . date() . "${LPR}run_bismark$N: WARNING: BAMFile $LCY$BAMFile$N is newer than fixedBAMFile $LGN$fixedBAMFile$N\n");
		push(@Ftorun, 5) if not grep (/^5$/, @Ftorun);
	}
	if (newer($BAMFile,$filteredBAMFile) == 1) {
		LOG($outLog, "\n" . date() . "${LPR}run_bismark$N: WARNING: BAMFile $LCY$BAMFile$N is newer than filteredBAMFile $LGN$filteredBAMFile$N\n");
		push(@Ftorun, 6) if not grep (/^6$/, @Ftorun);
	}
	if (newer($fixedBAMFile,$filteredBAMFile) == 1) {
		LOG($outLog, "\n" . date() . "${LPR}run_bismark$N: WARNING: BAMFile $LCY$fixedBAMFile$N is newer than filteredBAMFile $LGN$filteredBAMFile$N\n");
		push(@Ftorun, 6) if not grep (/^6$/, @Ftorun);
	}
	LOG($outLog, "\n$new1/$new2/$new3: Consider rerunning footLoop with$LPR -F " . join(",", @Ftorun) . "$N\n\n") if @Ftorun > 0;
	return("$outDir/$BAMFilename");
}

sub checkdate {
	my ($file) = @_;
	my ($date) = `date -r $file +\%s`;
	chomp($date);
	return($date);
}

sub newer {
	my ($file1, $file2) = @_;
	return(-1) if (checkdate($file1) < checkdate($file2));
	return( 1) if (checkdate($file1) > checkdate($file2));
	return 0;	
}

sub run_bismark_parallel {
	my ($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile, $SEQ, $outReadLog, $bismark_geneIndexDir, $bowtieOpt, $outLog) = @_;
	my $mapLog = "";
	$BAMFile =~ s/(.fq|.fastq|.fq.gz|.fastq.gz)_bismark_bt2/_bismark_bt2/;
	my $ext = "bam";
	my ($BAMFilename) = getFilename($BAMFile, "full");
	
	my $run_boolean = "\n\t${LGN}WAS NOT RUN$N:${YW} ";

	$run_boolean = "\n$YW\t ";

	############
	#  PARALEL #
	############

	my $outFolder = $outDir . "/.run_bismark/";
	system("mkdir -p $outFolder") if not -d $outFolder;
	

	if (defined $force{11} and -d $outFolder) {
		my @files = <$outFolder/*.part.gz>;
		if (@files != 0) {
			if (not defined $opt_0) {
				DIELOG($outLog, "DEBUG: Exited before removing *.part.gz in $LCY$outFolder/$N\n");
				system("/bin/rm $outFolder/*.part.gz");
			}
		}
	}

	##############################################
	# SplitFastq by -P reads each (default 1000) #
	##############################################
	my $footLoopSplitFastqScript = "$footLoop_script_folder/bin/footLoop_SplitFastq.pl";
	LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: Running ${LCY}$footLoopSplitFastqScript$N to split readFile into $LGN$parallel$N reads per file\n"); 
	my $splitFastq_resultFile = "$outFolder/.footLoop_SplitFastq_result.txt";
	my $splitFastq_cmd = $parallel ne -1 ? "$footLoopSplitFastqScript -i $readFile -o $outFolder -n $parallel > $splitFastq_resultFile" :
														"/bin/ln -s $readFile $outFolder/$readFilename.0.part.gz > $splitFastq_resultFile";
	print_cmd($splitFastq_cmd, $outLog);

	LOG($outReadLog, "footLoop.pl,run_bismark,$splitFastq_cmd\n","NA");

	my @fastqFiles = <$outFolder/*.part.gz>;
	my $totalfastqFiles = @fastqFiles;
	if (defined $force{3} or not -e $splitFastq_resultFile or $totalfastqFiles == 0) {
		if (not defined $opt_0) {
			system($splitFastq_cmd) == 0 or DIELOG($outLog, "footLoop.pl::run_bismark: Failed to run splitFastq_cmd: $LRD$!$N\n$LCY\n$splitFastq_cmd\n$N\n");
		}
		@fastqFiles = <$outFolder/*.part.gz>;
		$totalfastqFiles = @fastqFiles;
		LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: ${GN}SUCCESS$N: Fastq was split successfully into $LGN$totalfastqFiles$N .part.gz files in folder $LCY$outFolder/*.part.gz$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: ${GN}SUCCESS$N: Using previously run $footLoopSplitFastqScript $LGN$totalfastqFiles$N .part.gz files in folder $LCY$outFolder/*.part.gz$N\n");
	}

	if (-e $splitFastq_resultFile) {
		my $splitFastq_resultLog_head = -e $splitFastq_resultFile ? `head $splitFastq_resultFile` : "splitFastq_resultFile$N doesn't exist: $LCY$splitFastq_resultFile$N\n";
		my $splitFastq_resultLog_tail = -e $splitFastq_resultFile ? `tail $splitFastq_resultFile` : "splitFastq_resultFile$N doesn't exist: $LCY$splitFastq_resultFile$N\n";
		LOG($outLog, "\n$splitFastq_resultLog_tail$N\n","NA");
	}
	LOG($outLog, "==================================================\n");			

	##############################################
	# BISMARK ON EACH SPLIT FASTQ
	# IMPORTANT! --temp_dir FILENAME_bismark is important
	#   not doing this will cause all .fq temp files
	#   in the folder where the script is ran.. if multiple fq files with same name,
	#   then it'll be gg as they will overwrite each other
	#   e.g. sep.fastq.gz_C_to_T.fastq etc etc
	##############################################

	my $bismark_cmd = "bismark --non_directional -o FILENAME_bismark --temp_dir FILENAME_bismark --ambig_bam --ambiguous --unmapped $bowtieOpt $bismark_geneIndexDir FILENAME > FILENAME_bismark/.bismark_log 2>&1";

	LOG($outLog, "\n" . date() . "${LPR}bismark .part.gz .bam$N: Running ${LCY}bismark$N in parallel\n");
	print_cmd("$bismark_cmd", $outLog);

	my @bamFiles = <$outFolder/*_bismark/*bismark_bt2.bam>;
	my $totalbamFiles = @bamFiles;
	
	if (defined $force{4} or $totalfastqFiles > $totalbamFiles or $totalbamFiles == 0) {
		run_footLoop_sbatch($bismark_cmd, "bismark", "bam", \@fastqFiles, $max_slurm_job, $outLog);

		@bamFiles = <$outFolder/*_bismark/*bismark_bt2.bam>;

		# AMBIG
		if (defined $ambig) {
			for (my $i = 0; $i < @bamFiles; $i++) {
				my $ambigFile = $bamFiles[$i];
				$ambigFile =~ s/.bam$/.ambig.bam/;
				next if not -e $ambigFile;

				my $bamFile = $bamFiles[$i] . ".unmerged";
				system("/bin/mv $bamFiles[$i] $bamFile");
				my $fofnFile = "$bamFile.fofn";
				open (my $outfofn, ">", $fofnFile) or die "Failed to write to $fofnFile: $!\n";
				print $outfofn "$bamFile\n";
				print $outfofn "$ambigFile\n";
				close $outfofn;
				my $merge_bam_cmd = "samtools merge -f --reference $genomeFile $bamFiles[$i] -b $fofnFile";
				my $merge_bam_cmd_print = $merge_bam_cmd;
				print_cmd($merge_bam_cmd, $outLog);
				system($merge_bam_cmd) == 0 or DIELOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: ${LRD}FAILED$N: Failed to merge cmd: $!\nCMD=$LCY$merge_bam_cmd$N\n\n");
				LOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: Merging ambig into .bam files using samtools\n");
			}
		}




		$totalbamFiles = @bamFiles;
		LOG($outLog, "\n" . date() . "${LPR}bismark .part.gz .bam$N: ${GN}SUCCESS$N: Created $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$outFolder/*_bismark/$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}bismark .part.gz .bam$N: ${GN}SUCCESS$N: Using previously run $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$outFolder/*_bismark/$N\n");
	}
	LOG($outLog, "==================================================\n");			


	################
	# FIX BAM FILE #
	################
	LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: Running ${LCY}$fixBAMFile_script$N in parallel\n");

	my $fix_BAMFile_cmd = "$fixBAMFile_script -b FILENAME -g $seqFile -i $geneIndexFile -o FOLDER";
	print_cmd("$fix_BAMFile_cmd", $outLog);

	my @fixedFiles = <$outFolder/*_bismark/*.$ext.fixed.gz>;
	my $totalfixedFiles = @fixedFiles;

	if (defined $force{5} or $totalfixedFiles < $totalbamFiles) {
		run_footLoop_sbatch($fix_BAMFile_cmd, "fixBAM", "fixed.gz", \@bamFiles, $max_slurm_job, $outLog);
		@fixedFiles = <$outFolder/*_bismark/*.$ext.fixed.gz>;
		$totalfixedFiles = @fixedFiles;
		LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: ${GN}SUCCESS$N: Created $LGN$totalfixedFiles$N .fixed.gz fixed Bam Files in $LCY$outFolder/*_bismark/*.$ext\_fixBAM$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: ${GN}SUCCESS$N: Using previously run $LGN$totalfixedFiles$N .fixed.gz fixed Bam Files in $LCY$outFolder/*_bismark/*.$ext\_fixBAM$N\n");
	}
	LOG($outLog, "==================================================\n");			

	###################
	# FILTER BAM FILE #
	###################
	LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: Running ${LCY}$filterBAMFile_script$N in parallel\n");

	my $filter_BAMFile_cmd = "$filterBAMFile_script -r $opt_r -b FILENAME -g $seqFile -L $opt_L -O FOLDER -o FOLDER -i $geneIndexFile2";
	   $filter_BAMFile_cmd .= " -q $opt_q" if defined $opt_q;
	   $filter_BAMFile_cmd .= " -x $opt_x" if defined $opt_x;
	   $filter_BAMFile_cmd .= " -y $opt_y" if defined $opt_y;
	print_cmd("$filter_BAMFile_cmd", $outLog);

	my @filterBAMFile = <$outFolder/*_bismark/*.log.txt>; #log.txt coz there could be massive number of .filtered files which will slow down the script
	my $totalfilterBAMFile = @filterBAMFile;

	if (defined $force{6} or @fixedFiles < @bamFiles or $totalfilterBAMFile == 0) {
		run_footLoop_sbatch($filter_BAMFile_cmd, "filterBAM", "log.txt", \@bamFiles, $max_slurm_job, $outLog);
		@filterBAMFile = <$outFolder/*_bismark/*.log.txt>; 
		$totalfilterBAMFile = @filterBAMFile;
		LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: ${GN}SUCCESS$N: Created $LGN$totalfilterBAMFile$N filterBAM/*.log.txt files in $LCY$outFolder/*_bismark/*.$ext\_fixBAM or filterBAM/$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: ${GN}SUCCESS$N: Using previously run $LGN$totalfilterBAMFile$N filterBAM/*.log.txt files in $LCY$outFolder/*_bismark/*.$ext\_fixBAM or filterBAM/$N\n");
	}
	LOG($outLog, "==================================================\n");			
	
	##################
	# MERGE FIXED.GZ #
	##################
	my $output_merge_fixedFiles = "$outDir/$BAMFilename.fixed.gz";

	LOG($outLog, "\n" . date() . "${LPR}Merge *part/.fixed.gz .fixed.gz$N: Merging $LGN$totalfixedFiles$N .fixed.gz files into $LCY$output_merge_fixedFiles$N\n");

	my $merge_fixedFiles_cmd = "zcat -f $outFolder/*_bismark/*.$ext.fixed.gz | gzip > $output_merge_fixedFiles";
	print_cmd($merge_fixedFiles_cmd, $outLog);

	if (defined $force{7} or not -e $output_merge_fixedFiles) {
		if (not defined $opt_0) {
			system($merge_fixedFiles_cmd) == 0 or LOG($outLog, "\n" . date() . "${LPR}Merge *part/.fixed.gz .fixed.gz$N: Failed to merge_filtered_cmd: $!\n$LCY$merge_fixedFiles_cmd$N\n");
		}
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.fixed.gz .fixed.gz$N: ${GN}SUCCESS$N: Merged $LGN$totalfixedFiles$N .fixed.gz files into $LCY$output_merge_fixedFiles$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.fixed.gz .fixed.gz$N: ${GN}SUCCESS$N: Using previously made merged .fixed.gz $LCY$output_merge_fixedFiles$N\n");
	}
	LOG($outLog, "==================================================\n");			

	##########################################
	# Merge .filtered.gz into each gene file #
	##########################################
	my @strands = qw(Pos Neg);
	my $runtotal = scalar(keys %{$SEQ}) * @strands;
	my $runcount = 0;
	my $merge_success_created   = 0;
	my $merge_success_usingprev = 0;

	LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N: Merging .filtered.gz files into one file per gene\n");
	makedir("$outDir/origFiles") if not -d "$outDir/origFiles";
	my $merge_filtered_cmd_print = "zcat -f $outFolder/*_bismark/FILENAME | gzip -f > $outDir/origFiles/FILENAME";
	print_cmd($merge_filtered_cmd_print, $outLog);

	foreach my $gene (sort keys %{$SEQ}) {
		foreach my $strands (sort @strands) {
			LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N: Done $LGN$runcount/$runtotal$N\n") if $runcount % 100 == 0;
			$runcount ++;
			my $filteredFilename = "$gene\_$strands.filtered.gz";
			LOG($outLog, "$filteredFilename\n","NA");
			my $outputFile = "$outDir/origFiles/$filteredFilename";
			my $merge_filtered_cmd = "zcat -f $outFolder/*_bismark/$filteredFilename | gzip -f > $outputFile";
			if ($runcount == 1) {
				print_cmd($merge_filtered_cmd, $outLog);
			}
			LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N: $LCY$gene $strands$N: Merging .filtered file into $LCY$outputFile$N\n\n$LCY$merge_filtered_cmd$N\n","NA"); 
			if (defined $force{8} or not -e $outputFile) {
				if (not defined $opt_0) {
					my $merge_filtered_cmd_success = system($merge_filtered_cmd);
					if ($merge_filtered_cmd_success != 0) {
						LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N: ${LRD}FAILED$N: Failed to merge_filtered_cmd: $!\n$LCY$merge_filtered_cmd$N\n");
					}
					else {
						LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N: ${GN}SUCCESS$N: $LPR$runcount/$runtotal$N: Created .filter: $LCY$gene $LPR$strands$N: Created merged filtered $LCY$output_merge_fixedFiles$N\n","NA");
						$merge_success_created ++;
					}
				}
			}
			else {
				LOG($outLog, "\n" . date() . " ${GN}SUCCESS$N: $LPR$runcount/$runtotal$N: Using prev merged .filter: $LCY$gene $LPR$strands$N: Using previously made merged filtered $LCY$output_merge_fixedFiles$N\n","NA");
				$merge_success_usingprev ++;
			}
		}
	}
	if ($merge_success_created + $merge_success_usingprev == $runcount) {
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N ${GN}SUCCESS$N: Created $LGN$merge_success_created$N/$LGN$runtotal$N and using previously made $LGN$merge_success_usingprev$N/$LGN$runtotal$N .filtered.gz files (total=$LGN$runcount$N) $LCY$output_merge_fixedFiles$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.filtered.gz .filtered.gz$N ${LRD}FAIL$N: Created $LGN$merge_success_created$N/$LGN$runtotal$N and using previously made $LGN$merge_success_usingprev$N/$LGN$runtotal$N .filtered.gz files (total=$LGN$runcount$N) $LCY$output_merge_fixedFiles$N\n");
	}
	LOG($outLog, "==================================================\n");			
	
	#########################################
	# CREATE .fofn FILE FOR MERGE BAM FILES #
	#########################################
	my $fofnFile = "$outFolder/bamlist.fofn";
	my $fofn_cmd = "/bin/ls $outFolder/*_bismark/*bismark_bt2.$ext > $fofnFile";
	LOG($outLog, "\n" . date() . "${LPR}Create *part/.bam .fofn$N: Creating .fofn file\n");
	print_cmd($fofn_cmd, $outLog);
	system($fofn_cmd) == 0 or DIELOG($outLog, "Failed to write to $fofnFile: $!\n$LCY$fofn_cmd$N\n");
	my ($fofnSize) = `wc -l $fofnFile` =~ /^(\d+)\s*/;
	LOG($outLog, "\n" . date() . "${LPR}Create *part/.bam .fofn$N: ${GN}SUCCESS$N: Created fofn file of $LGN$fofnSize$N bam files ($LCY$fofnFile$N)!\n");
	LOG($outLog, "==================================================\n");			

	###################
	# MERGE BAM FILES #
	###################
	my @HEADER; my @REPORT;
	my $merge_bam_cmd = "samtools merge -f --reference $genomeFile $BAMFile -b $fofnFile";
	my $merge_bam_cmd_print = $merge_bam_cmd;
	LOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: Merging $LGN$totalbamFiles$N .bam files using samtools\n");
	print_cmd($merge_bam_cmd, $outLog);

	if (defined $force{9} or not -e $BAMFile) {
		if (not defined $opt_0) {
			system($merge_bam_cmd) == 0 or DIELOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: ${LRD}FAILED$N: Failed to merge cmd: $!\nCMD=$LCY$merge_bam_cmd$N\n\n");
		}
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: ${GN}SUCCESS$N: Merged $LGN$totalbamFiles$N .bam files into main .bam file $LCY$BAMFile$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}Merge *part/.bam .bam$N: ${GN}SUCCESS$N: Using previously merged $LGN$totalbamFiles$N .bam files into main .bam file $LCY$BAMFile$N\n");
	}
	LOG($outLog, "==================================================\n");			

	#########################
	# MERGE BISMARK REPORTS #
	#########################
	LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: Merging reports in outFolder=$LCY$outFolder*/*_SE_report.txt$N\n");
	print_cmd("#footLoop.pl::parse_bismark_report(\$part_bismark_report, \$outLog)", $outLog);
	my $report_not_found = 0;
	for (my $p = 0; $p < @bamFiles; $p++) {
		my $bamFiles = $bamFiles[$p];
		my ($partFilename) = $bamFiles =~ /^(.+).bam$/;
		my $part_bismark_report = $partFilename . "_SE_report.txt";
		if (not -e $part_bismark_report) {
			$report_not_found ++;
			LOG($outLog, "\nfootLoop.pl::run_bismark: $LGN$p/$totalbamFiles$N. can't find report file $LCY$part_bismark_report$N\n");
			next;
		}
		my ($header, $report) = parse_bismark_report($part_bismark_report, $outLog);
		my @header = @{$header};
		my @report = @{$report};
		@HEADER = @header if $p == 0;
		LOG($outLog, "\t#$p. $LCY$part_bismark_report$N:\n\t","NA") if $p < 5;
		for (my $q = 0; $q < @header; $q++) {
			# just next coz we're logging, not that important
			next if not defined $header[$q];
			next if $header[$q] ne $HEADER[$q];
			next if not defined $report[0];
			next if not defined $report[$q];
			$REPORT[$q] += $report[$q] if $header[$q] !~ /^perc_/;
			$REPORT[$q] += $report[0]*$report[$q]/100 if $header[$q] =~ /^perc_/; 
			LOG($outLog, "q=$q, $report[0]*$report[$q]/100 = $REPORT[$q]\n","NA") if $header[$q] =~ /^perc_/ and $p < 5;
		}
		LOG($outLog, "\n","NA") if $p < 5;
	}
	my $merged_report_to_log = "";
	my ($bismark_report) = $BAMFile =~ /^(.+).$ext/; $bismark_report .= "_footLoop_report.txt";
	my $bismark_report_exist = 0;
	if (-e $bismark_report) {
		$bismark_report_exist = 1;
	}
	open (my $outbismark_report, ">", $bismark_report) or DIELOG($outLog, "footLoop.pl::run_bismark: parallel: Failed to write to bismark report $LCY$bismark_report$N: $!\n");
	for (my $q = 0; $q < @HEADER; $q++) {
		$REPORT[$q] = int(10000*$REPORT[$q]/($REPORT[0])+0.5)/100 if $HEADER[$q] =~ /^perc_/; 
		print $outbismark_report "$HEADER[$q]\t$REPORT[$q]\n";
		$merged_report_to_log .= "$LCY$HEADER[$q]$N=$LGN$REPORT[$q]$N\n";
	}
	close $outbismark_report;

	$mapLog  = "footLoop.pl,map," . "header\t" . join("\t", @HEADER) . "\tfootLoop_outDir\tuuid\n";
   $mapLog .= "footLoop.pl,map," . "record\t" . join("\t", @REPORT) . "\t$outDir\t$uuid\n";
	if ($bismark_report_exist == 0) {
		LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: ${GN}SUCCESS$N: Merged $LGN$totalbamFiles$N bismark reports in $LCY$outFolder*/*_SE_report.txt$N, stats: $merged_report_to_log\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}Merge _SE_report.txt _footLoop_report.txt$N: ${GN}SUCCESS$N: Using previously merged $LGN$totalbamFiles$N bismark reports in $LCY$outFolder*/*_SE_report.txt$N, stats:\n\n$merged_report_to_log\n");
	}
	LOG($outLog, "\t#MAP_RESULT\n$LGN$mapLog$N\n");
	LOG($outLog, "==================================================\n");			


	# print to $outReadLog
	LOG($outReadLog, "footLoop.pl,BAMFile,$outDir/$BAMFilename\n","NA");
	LOG($outReadLog, $mapLog,"NA");

	return("$outDir/$BAMFilename");
}

sub bismark_check_fq_bam {
	my ($fastqFilesArr, $bamFilesArr, $outLog) = @_;
	if (not defined $fastqFilesArr or not defined $bamFilesArr) {
		my $totalFastq = 0 if not defined $fastqFilesArr;
		my $totalbamFiles = 0 if not defined $bamFilesArr;
		
		LOG($outLog, "\n" . date() . "ERROR: fastq_not_in_BAM:${LRD}ALL$N, totalFastq=$LGN$totalFastq$N, totalBAM=$LGN$totalbamFiles$N\n");
		return;
	}
	my @fastqFiles = @{$fastqFilesArr};
	my @bamFiles = @{$bamFilesArr};
	my %temp;
	for (my $i = 0; $i < @bamFiles; $i++) {
		$temp{BAM}{$bamFiles[$i]} = 1;
	}
	my $fastq_no_BAM = 0;
	for (my $i = 0; $i < @fastqFiles; $i++) {
		$temp{fastq}{$fastqFiles[$i]} = 1;
		my ($fastqFilesName) = getFilename($fastqFiles[$i], "full");
		my $fastqFilesBAMFile = $fastqFiles[$i] . "_bismark/$fastqFilesName\_bismark_bt2.bam";
		$fastq_no_BAM ++ if not defined $temp{BAM}{$fastqFilesBAMFile};
	}
	my $totalFastq = @fastqFiles;
	my $totalbamFiles = @bamFiles;
	LOG($outLog, "\n" . date() . "footLoop.pl::bismark_check_fq_bam: Fastq_without_BAM=$LGN$fastq_no_BAM$N,totalFastq=$LGN$totalFastq$N,totalBAM=$LGN$totalbamFiles$N\n");
}

sub run_footLoop_sbatch {
	my ($cmd, $suffix, $ext, $filesARRAY, $max_parallel_run, $outLog, $force_sbatch) = @_;

	# - suffix : $file\_$suffix.sbatch/sbout/folder
	# - ext    : determine donefile (<$file\_$suffix/*.$ext>)

	my $jobidhash;
	my $totalfiles = scalar(@{$filesARRAY});
	for (my $i = 0; $i < @{$filesARRAY}; $i++) {

		my $file = $filesARRAY->[$i];
		my ($folder, $filename) = getFilename($file, "folderfull");

		my $sbatchfile = $file . "\_$suffix.sbatch";
		my $sboutfile  = $file . "\_$suffix.sbout";
		my $donefile   = $file . "\_$suffix.done";

		my $cmdcopy = $cmd;
			$cmdcopy =~ s/FILENAME/$file/g;	

		if ($cmdcopy !~ /FOLDER/) {
			system("mkdir -p $file\_$suffix") if not -e "$file\_$suffix";
		}
		else {
			$cmdcopy =~ s/FOLDER/$folder/g;	
		}
		if ($i == 0) {
			print_cmd($cmdcopy, $outLog);
		}

		my $sbatchprint = "";
			$sbatchprint .= "#!/bin/bash -l\n";
			$sbatchprint .= "#SBATCH -n 2 -N 1 -p high --mem 8000 -t 999:99:99\n";
			$sbatchprint .= "#SBATCH --job-name \"$filename\_$suffix\"\n";
			$sbatchprint .= "#SBATCH --output \"$sboutfile\"\n\n";
			$sbatchprint .= "conda activate footLoop2\n" if defined $opt_8;
			$sbatchprint .= "$cmdcopy && echo \"Done!\" > $donefile\n\n";

		open (my $out, ">", $sbatchfile) or die "Can't write to $LCY$sbatchfile$N: $!\n";
		print $out $sbatchprint;
		close $out;

		if (-e $donefile) {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles run_footLoop_sbatch $suffix$N: sbatch $LCY$sbatchfile$N # ${LGN}DONE$N\n");
			next;
		}
		if (not -e $donefile) {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles run_footLoop_sbatch $suffix$N: sbatch: $LCY$sbatchfile$N\n");
		}
		
		if (defined $opt_0) { # Debug
			next;
		}

		if ($i != 0) {
			my $sleep = 0;
			while (1) {
				last if $i < $max_parallel_run;
				my ($job_left) = squeue_check($jobidhash);
				LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
				last if ($job_left < $max_parallel_run);
				$sleep ++;
				sleep 5;
			}
		}
		my ($jobid) = `sbatch $sbatchfile`;
		chomp($jobid);
		($jobid) = $jobid =~ /^Submi.+job (\d+)$/;
		next if not defined $jobid;
		$jobidhash->{$jobid} = 1;
		LOG($outLog, "$YW$i$N $LCY$filename$N $LGN$jobid$N\n");
	}
	my $sleep = 0;
	while (1) {
		my ($job_left) = squeue_check($jobidhash);
		LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
		last if $job_left == 0;
		$sleep ++;
		sleep 5;
	}
	LOG($outLog, "\n" . date() . "All have been run!\n\n");
	return(0);
}

sub squeue_check {
	my ($jobidhash, $outLog) = @_;
	my @squeue = `squeue`;
	my $squeuehash;
	foreach my $line (@squeue) {
		next if $line =~ /JOBID\s+PARTITION.+/;
		my ($jobid) = $line =~ /^\s*(\d+)\s+/;
		if (not defined $jobid) {
			LOG($outLog, "Can't parse jobid from line=$LCY$line$N\n");
			next; # just next so we don't kill the script...
		}
		next if not defined $jobidhash->{$jobid};
		$squeuehash->{$jobid} = 1;
	}
	foreach my $jobid (keys %{$jobidhash}) {
		next if defined $squeuehash->{$jobid};
		undef $jobidhash->{$jobid};
		delete $jobidhash->{$jobid};
	}
	my ($total) = scalar(keys %{$jobidhash});
	return ($total);
}
sub print_cmd {
	my ($cmd, $outLog) = @_;
	LOG($outBigCMD, "\n$cmd\n","NA");
	if ($cmd !~ /^#/) {
		$cmd =~ s/^/    /;
		$cmd =~ s/ \-/ \\\n      \-/g;
	}
	else {
		$cmd =~ s/^#/   # /;
		$cmd =~ s/ \-/ \\\n     # \-/g;
	}
	LOG($outLog, "$LGN\n$cmd\n$N\n");
	LOG($outBigCMD, "\n$cmd\n","NA");
}



sub parse_bismark_report {
	my ($bismark_report_file, $outLog, $is_footLoop_report) = @_;
	open (my $inz, "<", $bismark_report_file) or DIELOG($outLog, "footLoop.pl::parse_bismark_report: Failed to read from $LCY$bismark_report_file$N: $!\n");
	my @report; my @header;
	if (defined $is_footLoop_report) {
		while (my $line = <$inz>) {
			chomp($line);
			my ($header, $num) = split("\t", $line);
			push(@report, $num);
			push(@header, $header);
		}
	}
	else {
		while (my $line = <$inz>) {
			chomp($line);
			if ($line =~ /^Sequences analysed in total:/) {
				my ($num) = $line =~ /in total:[ \t]+(\d+)[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get total read from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "total_read");
			}
			if ($line =~ /^Mapping efficiency/) {
				my ($num) = $line =~ /efficiency:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Percent Mapping Efficiency from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "perc_mapped_read");
			}
			if ($line =~ /^Sequences which were discarded because genomic sequence could not be extracted/) {
				my ($num) = $line =~ /extracted:[ \t]+(\d+\.?\d*)[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Cannot Extract Chr from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "cannot_extract_chr");
			}
			if ($line =~ /^Total number of C's analysed:/) {
				my ($num) = $line =~ /Total number of C's analysed:[ \t]+(\d+\.?\d*)[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total C analysed from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "total_C");
			}
			if ($line =~ /^C methylated in CpG context:/) {
				my ($num) = $line =~ /C methylated in CpG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CpG from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "perc_meth_CpG");
			}
			if ($line =~ /^C methylated in CHG context:/) {
				my ($num) = $line =~ /C methylated in CHG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHG from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "perc_meth_CHG");
			}
			if ($line =~ /^C methylated in CHH context:/) {
				my ($num) = $line =~ /C methylated in CHH context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHH from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "perc_meth_CHH");
			}
			if ($line =~ /^C methylated in Unknown context:/) {
				my ($num) = $line =~ /C methylated in Unknown context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
				LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in Unknown from log $LCY$bismark_report_file$N!\n") if not defined $num;
				push(@report, $num);
				push(@header, "perc_meth_Unknown");
			}
		}
	}
	close $inz;
	my $mapLog  = "footLoop.pl,map," . "header\t" . join("\t", @header) . "\tfootLoop_outDir\tuuid\n";
	   $mapLog .= "footLoop.pl,map," . "record\t" . join("\t", @report) . "\t$outDir\t$uuid\n";
	return(\@header, \@report, $mapLog);
}


sub get_geneIndex_fasta {
	my ($geneIndexFile, $outDir, $logFile, $outLog) = @_;
	my $geneIndexName = getFilename($geneIndexFile);
	my $run_boolean = "\t${LGN}WAS NOT RUN$N:${YW} ";
	my $geneIndexFileNew = "$outDir/$geneIndexName\_$bufferL\_$bufferR\_bp.bed";
	LOG($outLog, "\n" . date() . "${LPR}bedtools_bed_change .bed _x_y.bed$N: Copying or transforming start and end coordinates (startBuffer $LGN$bufferL$N and endBuffer $LGN$bufferR$N) of ${LCY}geneIndexFile$N $LCY$geneIndexFile$N\n");
	if ($bufferL eq 0 and $bufferR eq 0) {
		my $buffer_cmd = "/bin/cp $geneIndexFile $geneIndexFileNew";
		print_cmd($buffer_cmd, $outLog);
		system($buffer_cmd) == 0 or DIELOG($outLog, "\tbedtools_bed_change: Failed because of $LGN$!$N cmd=$LCY$buffer_cmd$N\n\n");
		LOG($outLog, "\n" . date() . "${LPR}bedtools_bed_change .bed _x_y.bed$N: ${GN}SUCCESS$N: Copied ${LCY}geneIndexFile$N $LCY$geneIndexFileNew$N\n");
	}
	else {
		my $buffer_cmd = "$footLoop_script_folder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1";
		print_cmd($buffer_cmd, $outLog);
		system($buffer_cmd) == 0 or DIELOG($outLog, "\tbedtools_bed_change: Failed because of $LGN$!$N cmd=$LCY$buffer_cmd$N\n\n");
		LOG($outLog, "\n" . date() . "${LPR}bedtools_bed_change _x_y.bed$N: ${GN}SUCCESS$N: Created ${LCY}geneInexFile$N $LCY$geneIndexFileNew$N\n");
	}
	LOG($outLog, "==================================================\n");			

	return($geneIndexFileNew);
}

sub parse_geneIndexFile {
	my ($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog) = @_;
	my $geneIndexlinecount = 0;
	open (my $geneIndexInCheck, "<", $geneIndexFile) or DIELOG($outLog, "footLoop.pl::parse_geneIndexFile: Cannot read from ${LCY}geneIndexFile$N $geneIndexFile: $!\n");
	while (my $line = <$geneIndexInCheck>) {
		chomp($line); $geneIndexlinecount ++;
		next if $line =~ /^#/;
		my @arr = split("\t", $line);
		if (@arr < 6) {
			DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: ${LCY}geneIndexFile$N (-i $geneIndexFile) has to be a 6 column bed format! Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[1] !~ /^\d+$/) {
			DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: ${LCY}geneIndexFile$N (-i $geneIndexFile) column 1 isn't integer! ($arr[1]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[2] !~ /^\d+$/) {
			DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: ${LCY}geneIndexFile$N (-i $geneIndexFile) column 2 isn't integer! ($arr[2]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[4] !~ /^\-?\d+\.?\d*e?\-?\d*\.?\d*$/i) {
			DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: ${LCY}geneIndexFile$N (-i $geneIndexFile) column 4 isn't numeric ($arr[4]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[5] !~ /^[\+\-]$/) {
			DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: ${LCY}geneIndexFile$N (-i $geneIndexFile) column 6 isn't strand (- or +)! Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
	}
	close $geneIndexInCheck;

	$geneIndexFile = get_geneIndex_fasta($geneIndexFile, $outDir, $logFile, $outLog);
	LOG($outReadLog, "footLoop.pl,geneIndexFile,$geneIndexFile\n","NA");
	my $geneIndex; my $linecount = 0;
	open (my $geneIndexIn, "<", $geneIndexFile) or die "Cannot read from $geneIndexFile: $!\n";
	LOG($outLog, "\n" . date() . "${LPR}parse_geneIndexFile .bed .bed$N: Parsing genes from ${LCY}geneIndexFile$N $LCY$geneIndexFile$N\n");
	print_cmd("#footLoop.pl::parse_geneIndexFile", $outLog);
	while (my $line = <$geneIndexIn>) {
		chomp($line); $linecount ++;
		my ($chr, $beg, $end, $gene) = split("\t", $line);
		$gene = uc($gene);
		LOG($outLog, "\t$GR$linecount: gene=$LGN$gene$N,beg=$LGN$beg$N,end=$LGN$end$N,bufferLen=$LGN$bufferLen$N,length=$LGN" . ($end-$beg-$bufferLen) . "$N\n","NA");
		$geneIndex->{$gene} = $beg;
	}
	close $geneIndexIn;
	LOG($outLog, "\n" . date() . "${LPR}parse_geneIndexFile .bed .bed$N: ${GN}SUCCESS$N: Parsed $LGN$linecount$N gene coordinates from ${LCY}geneIndexFile$N $LCY$geneIndexFile$N\n");
	LOG($outLog, "==================================================\n");			

	my $geneIndexName = getFilename($geneIndexFile, "full");
	my $geneIndexFaDir = $outDir . "/.geneIndex/";
	makedir($geneIndexFaDir);
	my $geneIndexFaTemp = $outDir . "/.geneIndex/$geneIndexName.fa";

	LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: Getting fasta sequence using ${LCY}fastaFromBed$N\n");
	my $fastaFromBed_cmd = "fastaFromBed -fi $genomeFile -bed $geneIndexFile -fo $geneIndexFaTemp -nameOnly";
	print_cmd($fastaFromBed_cmd, $outLog);
	my $fastaFromBed_cmd_success = system("$fastaFromBed_cmd");
	if ($fastaFromBed_cmd_success == 0) {
		LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: ${GN}SUCCESS$N: Created fasta file $LCY$geneIndexFaTemp$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: ${LRD}FAIL$N: Failed to create fasta file $LCY$geneIndexFaTemp$N: $LGN$!$N\n");
	}
	LOG($outLog, "==================================================\n");			

	my ($geneIndexFaMD5, $temp, $geneIndexFaTempMD5File)  = getMD5($geneIndexFaTemp);
	LOG($outLog, "\tFatal error as MD5 of $geneIndexFaTemp is NA..(GENE=$geneIndexFaTemp, MD5=$geneIndexFaMD5)\n") and exit 1 if $geneIndexFaMD5 eq "NA";
	my $geneIndexFa        = $outDir . "/.geneIndex/$geneIndexFaMD5/$geneIndexName.fa";
	my $geneIndexFaMD5File = $outDir . "/.geneIndex/$geneIndexFaMD5/$geneIndexName.fa.md5";
	if (not -d $geneIndexFaMD5 or (-d $geneIndexFaMD5 and not -e $geneIndexFa)) {
		makedir("$outDir/.geneIndex/$geneIndexFaMD5");
		system("/bin/mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/") == 0 or LOG($outLog, "Failed to mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/: $!\n") and exit 1;
		system("/bin/mv $geneIndexFaTempMD5File $outDir/.geneIndex/$geneIndexFaMD5/") == 0 or LOG($outLog, "Failed to mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/: $!\n") and exit 1;
	}
	if (not defined $geneIndex or not defined $geneIndexFa) {
		DIELOG($outLog, "\nfootLoop.pl::parse_geneIndexFile: ERROR: geneIndexHash is not defined!\n");
	}
	my $SEQ = parse_fasta($geneIndexFa, $outLog, $minReadL);
	LOG($outReadLog, "footLoop.pl,-g,$genomeFile\n","NA");
	LOG($outReadLog, "footLoop.pl,seqFile,$geneIndexFa\n","NA");
	return ($SEQ, $geneIndex, $geneIndexFa, "$outDir/.geneIndex/$geneIndexFaMD5/", $geneIndexFile);
}

sub parse_fasta {
	my ($seqFile, $outLog, $minReadL) = @_;
	LOG($outLog, "\n" . date() . "${LPR}footLoop.pl::parse_fasta .fa .fa$N: Parsing reference sequences and infos from ${LCY}seqFile$N: $LCY$seqFile$N\n");
	open(my $SEQIN, "<", $seqFile) or DIELOG($outLog, "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!\n");
	my $fasta = new FAlite($SEQIN);
	my $linecount = 0;
	LOG($outLog, "\t${GR}From fasta file:$N\n","NA");
	while (my $entry = $fasta->nextEntry()) {
		$linecount ++;
		my $def = $entry->def; $def =~ s/^>//;
		if ($def =~ /^>.+::.+:\d+\-\d+$/) { #newer version of fastaFromBed
			($def) = $def =~ /^(.+)::.+:\d+\-\d+$/;
		}
		my $gene = uc($def);
		my $seqz = uc($entry->seq);
	   $SEQ->{$gene}{seq}       = [split("", $seqz)];
		$SEQ->{$gene}{geneL}     = scalar(@{$SEQ->{$gene}{seq}}) - $bufferLen;
		$SEQ->{$gene}{geneL}     = scalar(@{$SEQ->{$gene}{seq}}) if $SEQ->{$gene}{geneL} <= 0;
		$SEQ->{$gene}{minReadL}  = (defined $minReadL and $opt_L =~ /p$/i) ? int(0.5+$SEQ->{$gene}{geneL} * $minReadL / 100) : $minReadL;
		$SEQ->{$gene}{total}     = 0;
		$SEQ->{$gene}{badlength} = 0;
		$SEQ->{$gene}{lowq}      = 0;
		$SEQ->{$gene}{used}      = 0;
		$SEQ->{$gene}{pos}       = 0;
		$SEQ->{$gene}{neg}       = 0;
		$SEQ->{$gene}{orig}      = $def;
		LOG($outLog, "\t\t$linecount:gene=$LGN$gene$N,length=$LGN$SEQ->{$gene}{geneL}$N\n","NA");
	}
	close $SEQIN;
	LOG($outLog, "\n" . date() . "${LPR}footLoop.pl::parse_fasta .fa .fa$N: ${GN}SUCCESS$N: Parsed $LGN$linecount$N total sequences from fasta file\n");
	LOG($outLog, "==================================================\n");			

	return($SEQ);
}

sub uppercaseFasta {
	my ($fastaFile, $outLog) = @_;
	open (my $in, "<", $fastaFile) or LOG($outLog, "Failed to open $fastaFile: $!\n") and exit 1;
	open (my $out, ">", "$fastaFile.out") or LOG($outLog, "Failed to write to $fastaFile.out: $!\n") and exit 1;
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$seq =~ tr/a-z/A-Z/;
		print $out "$def\n$seq\n";
	}
	close $in;
	close $out;
	system("/bin/mv $fastaFile.out $fastaFile") == 0 or LOG($outLog, "Failed to /bin/mv $fastaFile.out $fastaFile: $!\n") and exit 1;
	return("$fastaFile");
}

sub check_sanity {
	my ($opts) = @_;

	my ($usageshort, $usage, $usagelong, $example) = getUsage();

	print $example and exit if defined $opt_e;
	my $checkopts = 0;
	if (defined $opts) {
		my %opts = %{$opts};
		foreach my $key (sort keys %opts) {
			$checkopts ++ if defined $opts{$key};
		}
	}
	die "\n$usageshort\n" if $checkopts == 0;
	die "\n$usage\n"      if defined $opt_h;
	die "\n$usagelong\n"  if defined $opt_H;

	my $errors = "";
	my $read0 = defined($opt_r) ? $opt_r : "FALSE";

	$errors .= "-r $LGN<sequencing_reads.fastq|fq|fastq.gz|fq.gz>$N: is not defined.\n" if not defined($opt_r);
	$errors .= "-r $LGN<sequencing_reads.fastq|fq|fastq.gz|fq.gz>$N: is defined ($opt_r) but file does not exist!\n" if defined $opt_r and not -e $opt_r;

	if (defined $opt_r and -e $opt_r) {
		$opt_n = $opt_r . "_footLoop" if not defined $opt_n;
		$opt_n = `realpath $opt_n`;chomp($opt_n);
		if (not -d $opt_n) {
			makedir($opt_n);
		}

		$errors .= "-L $LGN<min read length (<integer>:in bp; <integer>p:% amplicon>$N must be positive integer!\n" if defined($opt_L) and ($opt_L =~ /^0+\.?0*[p]?$/ or $opt_L !~ /^\d+[p]?$/);

		if (defined $opt_g and -e $opt_g and not defined $opt_i) {
			my (@faPaths) = split("/", $opt_g);
			my ($faFile) = $faPaths[@faPaths-1];
			my $prev_opt_g = `realpath $opt_g`;
			chomp($prev_opt_g);
	
			my $curr_opt_g = "$opt_n/$faFile";
			my $faiFile    = "$opt_n/$faFile.fai";
			$opt_i         = "$opt_n/$faFile.fai.bed";
	
			if (not -e $curr_opt_g) {
				my ($ln_success) = system("/bin/ln -s $prev_opt_g $opt_n/$faFile");
			}
			if (not -e $curr_opt_g) {
				print "footLoop.pl::check_sanity: can't create symlink of fasta file into output folder: $!\n$LCY$opt_g$N\n$LGN$opt_n$N\n\n";
				$opt_g = $prev_opt_g;
				$faiFile = $prev_opt_g . ".fai";
				$opt_i = $prev_opt_g . ".fai.bed";
			}
			else {
				$opt_g = $curr_opt_g;
				$faiFile = $curr_opt_g . ".fai";
				$opt_i = $curr_opt_g . ".fai.bed";
			}
			if (not -e $opt_i) {
				my $cmd = "samtools faidx $opt_g";
				my $cmdfai = system($cmd);
				$errors .= "-i not defined, tried to create $LGN$faiFile$N but can't\n" if $cmdfai ne 0;
				if (-e $faiFile) {
					fai_to_bed($faiFile, $opt_i);
					$errors .= "-i not defined, created faiFile $LCY$faiFile$N, and tried to create $LGN<$opt_i>$N but can't\n" if not -e $opt_i;
				}
			}
		}
		$errors .= "-g $LGN<ref_genome.fa>$N: is not defined.\n" if not defined($opt_g);
		$errors .= "-g $LGN<ref_genome.fa>$N: is defined ($YW$opt_g$N) but does not exist!\n" if defined $opt_g and not -e ($opt_g);
		$errors .= "-i $LGN<geneindex.bed>$N: is not defined.\n" if not defined($opt_i);
		$errors .= "-i <$LGN<geneindex.bed>$N: is defined ($YW$opt_i$N) but file does not exist!\n" if defined $opt_i and not -e ($opt_i);
	
		if (not -d "$footLoop_script_folder/.sortTMP/") {
			system("mkdir $footLoop_script_folder/.sortTMP/") == 0 or $errors .= "Failed to make directory $footLoop_script_folder/.sortTMP/: $!\n";
		}
	}
	
	if ($errors ne "") {
		print "$usageshort";
		print "\n${LRD}########## ${YW}FATAL ERROR${LRD} ##########\n\nREASON:$N\n";
		print $errors;
		print "\n\n${LRD}##################################$N\n\n" ;
		print "Use -Z if you're running invitro plasmid data!\n";
		die "\n";
	}
	print "Use -Z if you're running invitro plasmid data!\n";
}

sub fai_to_bed {
	my ($file, $outFile) = @_;
	$outFile = "$file.index.bed" if not defined $outFile;
	open (my $in1, "<", $file) or die "Failed to read from $LCY$file$N: $LRD$!$N\n\n";
	open (my $out1, ">", $outFile) or die "Failed to write to $outFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		my ($chr, $end) = split("\t", $line);
		die "Can't parse line!\n" if not defined $chr or not defined $end;
		my $beg = 0;
		print $out1 "$chr\t$beg\t$end\t$chr\t0\t+\n";
	}
	close ($in1);
	close $out1;
	print "$file\n$outFile\n\n";
}
sub record_options {
	my ($opts, $opts2, $outReadLog, $version, $outLog) = @_;
	my $optPrint = "$0";
	foreach my $opt (sort keys %{$opts}) {
		if (defined $opts->{$opt} and $opts->{$opt} eq "NONE") {
			$optPrint .= " -$opt" if defined $opts2->{$opt};
			LOG($outReadLog, "footLoop.pl,-$opt,\$TRUE\n","NA") if defined $opts2->{$opt};
			LOG($outReadLog, "footLoop.pl,-$opt,\$FALSE\n","NA") if not defined $opts2->{$opt};
		}
		elsif (defined $opts->{$opt} and $opts->{$opt} ne "NONE") {
			$optPrint .= " -$opt $opts->{$opt}";
			LOG($outReadLog, "footLoop.pl,-$opt,$opts->{$opt}\n","NA");
		}
	}
	my  $param = "
${YW}Initializing...$N

Version    : $version
Date       : $date
Run ID     : $uuid
Run script : $optPrint
	
${YW}Input Parameters$N
1. -r ${CY}Read/BAM$N  :$readFile
2. -g ${CY}Genome$N    :$genomeFile
3. -n ${CY}OutDir$N    :$outDir
4. -i ${CY}Index$N     :$geneIndexFile
5. -L ${CY}MinRdLen$N  :$minReadL
6. -q ${CY}MinMapQ$N   :$minMapQ
7. -Z ${CY}Bismark$N   :$bowtieOpt

";
	LOG($outLog, $param);
}

sub check_if_result_exist {
	my ($resFiles, $outLog) = @_;
	return if not defined $resFiles;
	my @resFiles = @{$resFiles};	
	return if @resFiles == 0;
	if (defined $opt_n and -d $opt_n and not defined $opt_F) {
		if (@{$resFiles} != 0) {
			my $resFileCheck = 0;
			foreach my $resFile (@{$resFiles}) {
				$resFileCheck = 1 if -e $resFile and -s $resFile > 0;
			}
			LOG($outLog, "-d $LGN$opt_n$N $resFiles->[0] already exist!\n") if $resFileCheck == 1;
			exit 0 if $resFileCheck == 1;
		}
		else {
			LOG($outLog, "End of run file does not exist!\n");
		}
	}
	return;
}


sub getUsage {
my $usageshort = "

-----------------
$YW $0 $version $N
-----------------

Usage: $YW$0$N -r <inputReads.fq.gz> -g <hg38.fa/plasmid.fa> [options..]

Required:
	-r input.fq (FASTQ file), accepts fastq/fastq.gz/fq/fq.gz

	-g genome.fa (FASTA file) e.g.$LGN
		>chr1
		ACTGGTGCA$N

Optional:
	-i geneIndex (bed6 or .fai from samtools faidx of -g <genome.fa>)
		#If not supplied, it'll create fasta index file from genome.fa and turn into bed6
		#bed6 example (tab separated):$LGN
		chr1\tstart\tend\tmygene1\t0\t+
		pCAG20\t0\t2390\tmygene2\t0\t+$N

	-n output directory [<input.fq.gz>_bismark]

Other options:
	-L [95p] minimum read length (add 'p' to indicate percent of fasta sequence length, e.g. seq length 1000bp, 95p means 950bp)
	-q [0] minimum map quality (use 0 since we're expecting a lot of conversions)
	-Z toggle to use non-stringent mapping 
		bowtie2 option --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 instead of
		default option --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
	-F toggle to redo bismark mapping even if a .BAM or .bam file is present in output_dir

	To prevent strange bismark error \"Chromosomal sequence could not be extracted\"
	-x [0] left  buffer in bp e.g. -x -10
	-y [0] right buffer in bp e.g. -y  10
	

Do $YW$0$N $LGN-h$N for longer explanation
Do $YW$0$N $LGN-e$N for example run

${LRD}IMPORTANT!!$N If you see a lot of 'Chromosomal sequence could not be extracted for..' try adding $YW-x -10 -y 10$N
-> If you still see then try adding bigger buffer, e.g. $YW-x -50 -y 50$N

-F:  (comma separated e.g. 2,3,4)
2: delete all part.gz
3: split Fastq into part.gz
4: bismark part.gz
5: filter bam ($fixBAMFile_script)
6: log bam ($filterBAMFile_script)
7: merge fixed.gz
8: merge filtered.gz
9: merge bam

";

my $usage = "

$LGN-Z$N: Make mapping parameters not stringent. 
Bismark uses bowtie2 mapping parameters:
- normal stringency (default) formula: --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
- lowest stringency (-Z)      formula: --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8

-H: Longer explanation about bismark parameters above.

(..continued; Use -H to display longer help!)

";

my $usagelong = "$usage

${YW}Required:$N
1) -r: Input Fastq file
2) -n: Output directoy name (will create if not exist)
3) -g: Fasta file of reference genome (e.g. hg19.fa)
4) -i: index bed file (bed6 with gene name on each)

${YW}Optional [default]:$N

-x: ${LGN}[0]$N Add this number from the$GN start$N of the index (strand taken into account)

-y: ${LGN}[0]$N Add this number from the$GN end$N of the index, for example:
\tcoordinate is:$GN chr1 200 500 SEQ 0 +$N and$GN -x -100 -y 50$N becomes: chr1 100 550 SEQ 0 +$N
\tcoordinate is:$GN chr1 200 500 SEQ 0 -$N and$GN -x -100 -y 50$N becomes: chr1 150 600 SEQ 0 -$N

-p: (not implemented) Run bismark script in parallel.

-c: consider Cs in CpG context (default: off; Don't include Cs in CpG context)

-L: ${LGN}[95p]$N minimum ${CY}read$N length in bp to be considered valid read
    Add \"p\" to make it percent of amplicon length
    e.g. 95p = 95% of amplicon length.

-F: Force re-create BAM file

-q: ${LGN}[0]$N q = map quality, where probability of mapped incorrectly = 10^(q/-10). 
Some quick -q references:
- 0  = 100  %
- 10 = 10   %
- 20 = 1    %
- 30 = 0.1  %
- 40 = 0.01 %

${LCY}PS: Using -q 0 when mapping barcoded reads onto KNOWN unique amplicon sequences will not produce incorrect result$N.


FYI, Bowtie2 mapping parameters used in bismark2:

--rdg <int1>,<int2>      Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets a penalty
                         of <int1> + N * <int2>. Default: 5, 3.

--rfg <int1>,<int2>      Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of length N gets
                         a penalty of <int1> + N * <int2>. Default: 5, 3.

--score_min <func>       Sets a function governing the minimum alignment score needed for an alignment to be considered
                         \"valid\" (i.e. good enough to report). This is a function of read length. For instance, specifying
                         L,0,-0.2 sets the minimum-score function f to f(x) = 0 + -0.2 * x, where x is the read length.
                         See also: setting function options at http://bowtie-bio.sourceforge.net/bowtie2. The default is
                         L,0,-0.2.
(bowtie2)
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)

Stringent         : --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
Not stringent (-Z): --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8
--rdg 5,3: read gap of length N penalty = 5+3N (stringent) or 2+N (not stringent)
--rfg 5,3: ref. gap of length N penalty = 5+3N (stringent) or 2+N (not stringent)
--score_min: minimum threshold of a read with length of L: -0.3L (stringent) or -0.8L (not stringent)

";

my $example = "

${GN}Example:$N

## Download hg19.fa.gz from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz, gunzip it, then put into footLoop folder.
## Download example.tar.gz and move everything in it into footLoop folder.
tar zxvf example.tar.gz
mv example/* footLoop/
cd footLoop
./footLoop.pl -r sample.fq.gz -n sample_MAP -g hg19.fa
./footPeak.pl -n sample_MAP -g hg19.fa -o sample_PEAK -t 55 -w 20 -l 50
./footClust.pl -n sample_PEAK
./footPeak_graph.pl -n sample_PEAK -r 1

";

	return($usageshort, $usage, $usagelong, $example);
}

