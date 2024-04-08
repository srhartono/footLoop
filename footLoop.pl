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

use warnings FATAL => 'all'; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars   qw($opt_v $opt_r $opt_g $opt_G $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f $opt_l $opt_e $opt_z $opt_9 $opt_o $opt_J $opt_0);
my @opts = qw(       $opt_r $opt_g $opt_i $opt_G $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_l $opt_e $opt_z $opt_9 $opt_o $opt_J $opt_0);
getopts("vr:g:G:i:n:L:x:y:q:HhZF:p:l:ez9o:J:0");

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

#my $version = "";
my %OPTS  = ('r' => $opt_r, 'g' => $opt_g, 'i' => $opt_i, 'n' => $opt_n, 'G' => $opt_G,
				 'L' => $opt_L, 'x' => $opt_x, 'y' => $opt_y, 'p' => $opt_p, 'J' => $opt_J,
				 'q' => $opt_q, 'F' => 'NONE', 'Z' => 'NONE', 'o' => $opt_o, '0' => $opt_0,
				 'p' => 'NONE', 'q' => $opt_q, 'l' => $opt_l);
my %OPTS2 = ('p' => $opt_p, 'Z' => $opt_Z, 'F' => $opt_F);

my $footLoop_script_folder = dirname(dirname abs_path $0) . "/footLoop/";
my $version = "unk";
my $md5script = "unk";
check_sanity(\%OPTS, \%OPTS2);
#my ($footLoop_script_folder, $version, $md5script) = check_software();
($footLoop_script_folder, $version, $md5script) = check_software();

print "version $YW$version$N\n";
if (not defined $opt_p) {
	$opt_p = -1;
}

###################
# 1. Define Input #
###################
my $snakemake;
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
#if (defined $parallel and $parallel ne -1) {# and $parallel < 20) {
#	print "Are you sure you'll divide read sequence by $parallel reads each? (press any key to continue)\n";
#	<STDIN>;
#}


my %force;
if (defined $opt_F) {
	my @forces = split(",", $opt_F);
	foreach my $force (@forces) {
		$force{$force} = 1;
		#print "force{$force}\n";
	}
	if (defined $force{1}) {print "1: rerun all fasta preprocessin\n";}
	if (defined $force{2}) {print "2: delete all part.gz\n";}
	if (defined $force{3}) {print "3: split Fastq into part.gz\n";}
	if (defined $force{4}) {print "4: bismark part.gz\n";}
	if (defined $force{5}) {print "5: filter bam (footLoop_2_fixBAMFile.pl)\n";}
	if (defined $force{6}) {print "6: log bam (footLoop_3_filterBAMFile.pl)\n";}
	if (defined $force{7}) {print "7: merge fixed.gz\n";}
	if (defined $force{8}) {print "8: merge filtered.gz\n";}
	if (defined $force{9}) {print "9: merge bam\n";}
}

#my ($label)    = defined $opt_l ? $opt_l : $outDir;#rigDir;# =~ /PCB[\-_0]*\d+/i ? $filteredDir =~ /(PCB[\-_0]*\d+)/i : $outDir;
#if (defined $opt_l) {
#	$label = $opt_l;
#}

#die "label=$label, uuid=$LCY$uuid$N\n";

# Make directory
makedir($outDir);

# make .LABEL
#open (my $outLabel, ">", "$outDir/.LABEL") or die "Failed to write to $outDir/.LABEL: $!\n";
#print $outLabel "$label\n";
#close $outLabel;
my $STEP = 0;

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
#($STEP) = LOGSTEP($outLog, "START", $STEP, 1, "Preprocessing gene index file and getting sequences\n");

# Add buffers to geneIndexFile ($footLoop_script_folder/lib/bedtools_bed_change.pl) and get their sequences using fastaFromBed, then get their info from fasta=$seqFile

my ($SEQ, $geneIndexHash, $seqFile, $bismark_geneIndexDir, $geneIndexFile2) = parse_geneIndexFile($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog);
#seqFile is fasta file that's inside the $outDir/.geneIndex/<MD5> folder, e.g. $outDir/.geneIndex/<MD5>/seqFile.fa
#	return ($SEQ, $geneIndex, $geneIndexFa, "$outDir/.geneIndex/$geneIndexFaMD5/");
#my $bismark_genome_preparation_folder = defined $opt_G ? $opt_G : "";

#LOGSTEP($outLog);

#####################################
# 2. Run Bismark Genome Preparation #
#####################################
#($STEP) = LOGSTEP($outLog, "START", $STEP, 1, "Creating bismark index\n");#$N $bowtieOpt $bismark_geneIndexDir $readFile\n");

# Make Bismark Index
bismark_genome_preparation($seqFile, $bismark_geneIndexDir, $bowtieOpt, $outLog);

#LOGSTEP($outLog);

##################
# 3. Run Bismark #
##################
#($STEP) = LOGSTEP($outLog, "START", $STEP, 1, "Running bismark\n");#$N $bowtieOpt $bismark_geneIndexDir $readFile\n");

# Run Bismark
($BAMFile) = run_bismark($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile2, $SEQ, $outReadLog, $outLog);
my ($BAMFileName) = getFilename($BAMFile);


#LOGSTEP($outLog);

LOG($outLog, "\n" . date() . "footLoop ran successfuly\n");
exit 0;
#DIELOG($outLog, "DEBUG: Exit before 3. splitting_BAMFile\n");
################################
# 3. Fixing BAM File #
################################

if (defined $opt_9) {
	die "-9 Defined so stopping here!\n";
}
($STEP) = LOGSTEP($outLog, "START", $STEP, 1, "Fix BAM File\n");#$N $bowtieOpt $bismark_geneIndexDir $readFile\n");

# Do footLoop_2_fixBAMFile.pl
# - Determine strand of read based on # of conversion
# - Determine bad regions in read (indels) which wont be used in peak calling
($BAMFile, $filteredDir) = split_BAMFile($BAMFile, $seqFile, $outReadLog, $outLog); #becomes .fixed
LOG($outReadLog, "footLoop.pl,filteredDir,$filteredDir\n");

LOGSTEP($outLog);

DIELOG($outLog, "DEBUG: Exit before 4. parse_BAMFile\n");
#	DIELOG($outLog, "DEBUG GOOD\n");


################################
# 4. Parse and Filter BAM File #
################################
($STEP) = LOGSTEP($outLog, "START", $STEP, 1, "Parse and Filter BAM File\n");#$N $bowtieOpt $bismark_geneIndexDir $readFile\n");

parse_BAMFile($BAMFile, $seqFile, $filteredDir, $outReadLog, $outLog);

LOGSTEP($outLog);
DIELOG($outLog, "GOOD\n");
#print_R_heatmap($SEQ);

###############
# SUBROUTINES #
###############

sub split_BAMFile {
	my ($BAMFile, $seqFile, $outReadLog, $outLog) = @_;
	LOG($outLog, "\n\ta. Fixing BAM file $CY$BAMFile$N with $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl\n");
	my ($BAMMD5) = getMD5($BAMFile);
	
	my $filteredDir = "$outDir/.BAMFile_$BAMMD5/";
	check_if_result_exist(["$filteredDir/.GOOD"], $outLog);
	system("mkdir -p $filteredDir") if not -d $filteredDir;
	my $checkBAM = 1;
	my ($BAMFileName) = getFilename($BAMFile, "full");
	my $BAMFileGZ = "$filteredDir/$BAMFileName.fixed.gz";
	$checkBAM = 0 if not -e "$filteredDir/$BAMFileName.fixed" and not -e "$filteredDir/$BAMFileName.fixed.gz";
	if (-e "$filteredDir/$BAMFileName.fixed.gz") {
		my ($BAMLineCount2) = linecount("$filteredDir/$BAMFileName.fixed.gz");
		my ($BAMLineCount1) = linecount($BAMFile);
		$checkBAM = $BAMLineCount1 - 10 > $BAMLineCount2 ? 0 : 2;
		if ($checkBAM eq 0) {
			LOG($outLog, "\n" . date() . "\tfootLoop.pl subroutine split_BAMFile::$LGN SUCCESS!!$N fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N exists (MD5=$LGN$BAMMD5$N) and total row $LGN($BAMLineCount2)$N >= total BAMFile row $LGN($BAMLineCount1 - 500)$N ($LCY$BAMFile$N)!\n");
		}
	}
	if (-e "$filteredDir/$BAMFileName.fixed" and $checkBAM == 0) {
		my ($BAMLineCount2) = linecount("$filteredDir/$BAMFileName.fixed");
		my ($BAMLineCount1) = linecount($BAMFile);
		$checkBAM = $BAMLineCount1 - 10 > $BAMLineCount2 ? 0 : 1;
		if ($checkBAM eq 0) {
			LOG($outLog, "\n" . date() . "\tfootLoop.pl subroutine split_BAMFile::$LGN SUCCESS!!$N fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N exists (MD5=$LGN$BAMMD5$N) and total row $LGN($BAMLineCount2)$N >= total BAMFile row $LGN($BAMLineCount1 - 500)$N ($LCY$BAMFile$N)!\n");
		}
	}

	my $fixBAMFile_cmd = "$footLoop_script_folder/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir";
	if ($checkBAM == 0) {
		LOG($outLog, "\tfootLoop.pl subroutine split_BAMFile:: fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N or .gz does not exist!\n");
		LOG($outLog, "\n$LCY$fixBAMFile_cmd$N\n\n");
		#DIELOG($outLog, "DEBUG Exited before running fixBAMFile_cmd\n");
		system($fixBAMFile_cmd) == 0 or DIELOG($outLog, "\n\n" . date() . "Failed to run footLoop_2_fixBAMFile.pl: $LCY$!$N\n\n$LCY$fixBAMFile_cmd$N\n\n");
		#system("$footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir") == 0 or LOG($outLog, "Failed to run $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir: $!\n") and exit 1;
		LOG($outReadLog, "footLoop.pl,split_BAMFile,$footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir\n","NA");
		if (not -e "$filteredDir/$BAMFileName.fixed.gz") {
			LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed");
			system("gzip $filteredDir/$BAMFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $filteredDir/$BAMFileName.fixed: $!\n");
		}
		else {
			LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed\n\t${LGN}Already exist! So not overwriting!\n");
		}
		$checkBAM = 1;
	}
	else {
		LOG($outLog, "\t${LGN}WAS NOT RUN$N: ${YW}::: $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -s $seqFile -o $filteredDir :::$N\n");
		LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed\n\t${LGN}Already exist! So not overwriting!\n");
	}

	# re-md5 BAMfile gz
	LOG($outLog, "\t${YW}$md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5$N\n");
	if ($md5script eq "md5sum") {
		system("$md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	if ($md5script eq "md5") {
		my ($md5res) = `$md5script $BAMFileGZ` =~ /^.+\= (\w+)$/; die "Failed to $md5script $BAMFileGZ: $!\n" if not defined $md5res;
		system("echo '$md5res $BAMFileGZ' > $filteredDir/.$BAMFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	LOG($outReadLog, "footLoop.pl,BAMFixedFile,$BAMFileGZ\n","NA");
	LOG($outReadLog, "footLoop.pl,BAMFixedFileMD5,$BAMMD5\n","NA");
	return($BAMFile, $filteredDir);
}

sub parse_BAMFile {
	my ($BAMFile, $seqFile, $filteredDir, $outReadLog, $outLog) = @_;
	my ($BAMFileName) = getFilename($BAMFile);

	LOG($outLog, "\ta. Parsing BAM file $CY$BAMFile$N and getting only high quality reads\nfilename=$BAMFileName\n\n");
	open(my $notused, ">", "$outDir/.$readFilename.notused") or (LOG($outLog, "Cannot open $outDir/.$readFilename.notused: $!\n") and exit 1);
	my $BAM;
	#open($BAM, $BAMFile) or (LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $BAMFile: $!") and exit 1 if $BAMFile =~ /.sam$/);
	open($BAM, "samtools view $BAMFile|") or (LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $BAMFile: $!") and exit 1) if $BAMFile =~ /.bam$/;
	
	## Some stats
	my $linecount   = 0;
	my $BAMStats = makehash(['total','used','diffgene','lowq','badlength']);
	
	while(my $line = <$BAM>) {
		#if ($linecount == 2) {DIELOG($outLog, "DEBUG GOOD\n");}
		$linecount ++;
		chomp($line);
		
		#####################
		# 1. Parse BAM line #
		#####################
		my @arr = split("\t", $line);
		LOG($outLog, "\t$YW$BAMFile$N: Done $linecount\n") if $linecount % 5000 == 0;
	
		# a. Total column must be 14 or skipped (e.g. BAM header)
		if (@arr < 14) {
			LOG($outLog, "$LGN$linecount$N: BAM header detected (#column < 14). LINE:\n\t$line\n");
			next;
		}
		
		my ($eval, $evalPrint) = myeval(\@arr);
		my ($read, $readStrand, $gene, $readPos, $mapQ, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $tags, $junk4, $readMethCall) = @arr;
		my $seqLen = length($seqs);
		$gene = uc($gene);
		check_BAM_field($BAMFile, $outLog, @arr);
		LOG($outLog, "\t$BAMFile: length of gene $gene undef\n") and exit 1 if not defined($SEQ->{$gene}{geneL});

		$BAMStats->{total} ++;
		$BAMStats->{gene}{$gene} ++;
	
		# Filter out reads which gene isn't in genomeFile and indexFile
		if (not defined $SEQ->{$gene}) {
			LOG($outLog, "$LGN$linecount$N: read $LCY$read$N is mapped to gene $LGN$gene$N is not in geneIndex! LINE:\n\t$line\n") if not defined $SEQ->{$gene};
			next;
		}
		
		$SEQ->{$gene}{total} ++;
		my ($readname) = length($read) >= 20 ? $read =~ /(.{20})$/ : $read; $readname = "..." . $readname  if length($read) >= 20; $readname = $read if not defined $readname;
		LOG($outLog, "\tExample at read $BAMStats->{total}: name=$CY$readname$N\tstrand=$CY$readStrand$N\tchr/gene=$CY$gene$N\tpos=$CY$readPos$N\tmapQ=$CY$mapQ$N\n\n") if $BAMStats->{total} == 1;
		LOG($outLog, "\tDone $GN$BAMStats->{total}$N\n") if $BAMStats->{total} % 2000 == 0;


		# a. Filter out reads with mapping quality=$mapQ less than threshold quality=$opt_q
		if($mapQ < $minMapQ)	{
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'lowq');
			print $notused "\t$CY$read$N quality ($CY$mapQ$N) is less than $CY$opt_q$N\n";
			$BAMStats->{lowq} ++;
			$SEQ->{$gene}{lowq} ++;
			next;
		}
		
		# b. Filter out reads with read length=$seqLen
		elsif (parse_cigar($cigar, "len") < $SEQ->{$gene}{minReadL}) {
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'badlength');
			my $cigarLen = parse_cigar($cigar, "len");
			print $notused "\t$CY$read$N length of seq (cigarLen = $LGN$cigarLen$N, seqLen = $CY$seqLen$N) is less than $SEQ->{$gene}{minReadL} bp (length of original sequence is ($CY" . $SEQ->{$gene}{geneL} . "$N)!\n";
		}

		# c. Filter out reads with more than 5% indel
		elsif (count_indel($cigar) > 5) {
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'indel');
			my ($indelperc) = count_indel($cigar);
			
			print $notused "\t$CY$read$N has high indel ($indelperc)!\n";
		}

		# c. If read strand is 0 or 16 store it
		elsif ($readStrand == 0 || $readStrand == 16)	{
			$SEQ->{$gene}{read}{$read} = $readStrand;
		}
		
		# d. Otherwise put it into notused
		else {
			print $notused "\t$CY$read$N isn't used for some reason\n";
		}
	}
	LOG($outLog, "\n\tb.Logging BAM file $CY$BAMFile$N\n");
	log_BAMFile($SEQ, $BAMFile, $BAMStats, $filteredDir, $outReadLog, $outLog);
}

sub count_indel {
	my ($cigar) = @_;
	my ($nums, $alps, $lenz) = parse_cigar($cigar);
	my ($ins, $del, $mat, $oth, $tot) = (0,0,0,0,0);
	for (my $i = 0; $i < @{$nums}; $i++) {
		my $alp = $alps->[$i];
		my $num = $nums->[$i];
		$mat += $num if $alp eq "M";
		$ins += $num if $alp eq "I";
		$del += $num if $alp eq "D";
		$oth += $num if $alp !~ /^(M|I|D)$/;
		$tot += $num;
	}
	my $matperc = $tot == 0 ? 0 : int($mat/$tot*1000+0.5)/10;
	my $insperc = $tot == 0 ? 0 : int($ins/$tot*1000+0.5)/10;
	my $delperc = $tot == 0 ? 0 : int($del/$tot*1000+0.5)/10;
	my $othperc = $tot == 0 ? 0 : int($oth/$tot*1000+0.5)/10;
	return $insperc+$delperc;
}

sub check_BAM_field {
	my ($BAMFile, $outLog, @arr) = @_;
	my $check = 0;
	for (my $i = 0; $i < @arr; $i++) {
		if ( not defined $arr[$i]) {
			LOG($outLog, "\n" . date() . "footLoop.pl::parse_BAMFile::check_BAM_field $BAMFile: column $i undefined\n");
			$check = 1;
		}
	}
	DIELOG($outLog, "footLoop.pl::parse_BAMFile::check_BAM_field: $BAMFile: Something is wrong with the bam file!\n") if $check != 0;
}

sub log_BAMFile {
	my ($SEQ, $BAMFile, $BAMStats, $filteredDir, $outReadLog, $outLog) = @_;
	foreach my $genez (sort keys %{$SEQ}) {
		my $outTXTFilePos  = "$filteredDir/$genez\_Pos.filtered";
		my $outTXTFileNeg  = "$filteredDir/$genez\_Neg.filtered";
		#my $outTXTFileUnk  = "$filteredDir/$genez\_Unk.filtered";
		open ($SEQ->{$genez}{outTXTPos}, ">", $outTXTFilePos);
		open ($SEQ->{$genez}{outTXTNeg}, ">", $outTXTFileNeg);
		#open ($SEQ->{$genez}{outTXTUnk}, ">", $outTXTFileUnk);
	}
	my ($BAMFileName) = getFilename($BAMFile);
	
	my $skipped = 0; my ($passedFilterP, $passedFilterN) = (0,0);
	open (my $inBAMFix, "zcat $filteredDir/$BAMFileName.fixed.gz|") or LOG($outLog, "Failed to open $filteredDir/$BAMFileName.fixed.gz: $!\n") and exit 1;
	#open (my $inBAMFix, "zcat < $filteredDir/$BAMFileName.fixed.gz|") or LOG($outLog, "Failed to open $filteredDir/$BAMFileName.fixed.gz: $!\n") and exit 1;
	while (my $line = <$inBAMFix>) {
		chomp($line);
		my ($read, $type, $oldStrand, $strand, $genez, $pos, $info) = split("\t", $line);
		if (not defined $SEQ->{$genez} or not defined $info) {
			LOG($outLog, "\tfootLoop.pl: ERROR in $LCY$filteredDir/$BAMFileName.fixed.gz$N: gene=$genez but \$SEQ->{\$genez} is not defined!line=\n$line\n\n") if not defined $SEQ->{$genez};
			DIELOG($outLog, "\tfootLoop.pl: ERROR in $LCY$filteredDir/$BAMFileName.fixed.gz$N: gene=$genez and seq genez is $SEQ->{$genez} but info is not defined!line=\n$line\n\n") if defined $SEQ->{$genez} and not defined $info;
			$skipped ++;
			next;
		}
		if (not defined $SEQ->{$genez}{read}{$read}) {
			next;
		}
		#my ($CT0, $CC0, $GA0, $GG0, $CT1, $CC1, $GA1, $GG1) = split(",", $info);
		#if ($type eq "6_BOTH" or $type eq "3_NONE") {
		#	print {$SEQ->{$genez}{outTXTUnk}} "$read\tBP\t$pos\n" if $strand eq 0;
		#	print {$SEQ->{$genez}{outTXTUnk}} "$read\tBN\t$pos\n" if $strand eq 16;
		#	$SEQ->{$genez}{unkpos} ++ if $strand == 0;
		#	$SEQ->{$genez}{unkneg} ++ if $strand == 16;
		#	$passedFilterP ++ if $strand == 0;
		#	$passedFilterN ++ if $strand == 16;
		#}
		#else {
			print {$SEQ->{$genez}{outTXTNeg}} "$read\tFN\t$pos\n" if $strand eq 16;
			print {$SEQ->{$genez}{outTXTPos}} "$read\tFP\t$pos\n" if $strand eq 0;
			$SEQ->{$genez}{pos} ++ if $strand == 0;
			$SEQ->{$genez}{neg} ++ if $strand == 16;
			$passedFilterP ++ if $strand == 0;
			$passedFilterN ++ if $strand == 16;
		#}
		$BAMStats->{used} ++;
		$SEQ->{$genez}{used} ++;
		$SEQ->{$genez}{posneg} ++ if $strand eq 16 and $oldStrand eq 0;
		$SEQ->{$genez}{negpos} ++ if $strand eq 0 and $oldStrand eq 16;
	}
	close $inBAMFix;
	#LOG($outLog, "\t$LCY$filteredDir/$BAMFileName.fixed$N: skipped = $LGN$skipped$N\n");

	#LOG($outReadLog, "footLoop.pl,read_passed_filter,header\ttotal\tpositive\tnegative\n");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter,record\t$BAMStats->{total}\t$passedFilterP\t$passedFilterN\n");

	LOG($outLog, "
Reads that passed filters:
Positive: $passedFilterP
Negative: $passedFilterN
Total   : $BAMStats->{total};
	
Per Gene:
	","NA");
	
	my $zero = "";
	foreach my $gene (sort keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		foreach my $key (@key) {
			$SEQ->{$gene}{$key} = 0 if not defined $SEQ->{$gene}{$key};
		}
	}
	foreach my $gene (sort {$SEQ->{$b}{total} <=> $SEQ->{$a}{total}} keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		my $outTXTFilePos  = "$filteredDir/$gene\_Pos.filtered"; system("/bin/rm $outTXTFilePos") if -e $outTXTFilePos and -s $outTXTFilePos == 0;
		my $outTXTFileNeg  = "$filteredDir/$gene\_Neg.filtered"; system("/bin/rm $outTXTFileNeg") if -e $outTXTFileNeg and -s $outTXTFileNeg == 0;
		#my $outTXTFileUnk  = "$filteredDir/$gene\_Unk.filtered"; system("/bin/rm $outTXTFileUnk") if -e $outTXTFileUnk and -s $outTXTFileUnk == 0;
	
		$zero .= "$gene ($SEQ->{$gene}{total})\n" and next if $SEQ->{$gene}{total} <= 10;
		$SEQ->{$gene}{total} = 0 if not defined $SEQ->{$gene}{total};
		$SEQ->{$gene}{orig} = 0 if not defined $SEQ->{$gene}{orig};
		$SEQ->{$gene}{negpos} = 0 if not defined $SEQ->{$gene}{negpos};
		$SEQ->{$gene}{posneg} = 0 if not defined $SEQ->{$gene}{posneg};
		$SEQ->{$gene}{unkpos} = 0 if not defined $SEQ->{$gene}{unkpos};
		$SEQ->{$gene}{unkneg} = 0 if not defined $SEQ->{$gene}{unkneg};
		$SEQ->{$gene}{used} = 0 if not defined $SEQ->{$gene}{used};
		$SEQ->{$gene}{total} = 0 if not defined $SEQ->{$gene}{total};
		$SEQ->{$gene}{badlength} = 0 if not defined $SEQ->{$gene}{badlength};
		$SEQ->{$gene}{indel} = 0 if not defined $SEQ->{$gene}{indel};
		$SEQ->{$gene}{lowq} = 0 if not defined $SEQ->{$gene}{lowq};
		$SEQ->{$gene}{pos} = 0 if not defined $SEQ->{$gene}{pos};
		$SEQ->{$gene}{neg} = 0 if not defined $SEQ->{$gene}{neg};
		my $gene2 = $SEQ->{$gene}{orig};
		my $text = "
- $gene (original name = $gene2):
Positive    = $SEQ->{$gene}{pos} (neg->pos = $SEQ->{$gene}{negpos})
Negative    = $SEQ->{$gene}{neg} (pos->neg = $SEQ->{$gene}{posneg})
UnkPos      = $SEQ->{$gene}{unkpos}
UnkNeg      = $SEQ->{$gene}{unkneg}
Used        = $SEQ->{$gene}{used}
Total       = $SEQ->{$gene}{total}
Too Short   = $SEQ->{$gene}{badlength}
High Indel  = $SEQ->{$gene}{indel}
Low Quality = $SEQ->{$gene}{lowq}
";
	#LOG($outLog, $text,"NA");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,header\tBAMFilename\tgene\ttotal\tused\tpositive\tnegative\tunkpos\tunkneg\ttooshort\tlowqual\n");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,record\t$BAMFileName\t$gene\t$SEQ->{$gene}{total}\t$SEQ->{$gene}{used}\t$SEQ->{$gene}{pos}\t$SEQ->{$gene}{neg}\t$SEQ->{$gene}{unkpos}\t$SEQ->{$gene}{unkneg}\t$SEQ->{$gene}{badlength}\t$SEQ->{$gene}{lowq}\n\n");
	#LOG($outReadLog, "footLoop.pl,gene_skipped_low_read,$BAMFileName,$zero\n","NA");
	}
	
	$zero = $zero eq "" ? "(None)\n" : "\n$zero\n";

	#LOG($outLog, "\n(Genes that have <= 10 reads:) $LGN$zero$N\n","NA");
	LOG($outLog, "\n" . date() . " ${GN}SUCCESS$N: Total=$BAMStats->{total}, used=$BAMStats->{used}, Low Map Quality=$BAMStats->{lowq}, Too short=$BAMStats->{badlength}\n");
	#LOG($outLog, "Output: ${LGN}$filteredDir$N\n\n");

}
sub filter_BAM_read {
	my ($seq, $count, $readname, $gene, $reason) = @_;
	$count->{$reason} ++;
	$seq->{$gene}{$reason} ++;
	return($seq, $count);
}

sub get_lotsOfC {
	my ($seq) = @_;
	my @seq = $seq =~ /ARRAY/ ? @{$seq} : split("", $seq);
	my $len = @seq;
	my $data;
	my %nuc;
	foreach my $nuc (@seq[0..@seq-1]) {
		$nuc = uc($nuc);
		$nuc{$nuc} ++;
	}
	foreach my $nuc (sort keys %nuc) {
		$nuc = uc($nuc);
		next if $nuc{$nuc} < 10;
		next if $nuc !~ /^(C|G)$/i;
		while ($seq =~ /${nuc}{6,$len}/ig) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my ($beg, $end) = (length($prev), length($prev) + length($curr));
			$data .= ",$nuc;$beg;$end";
			my ($prev0) = length($prev) > 2 ? $prev =~ /(..)$/ : "NOPREV";
			my ($next0) = length($next) > 2 ? $next =~ /^(..)/ : "NONEXT";
		}
	}
	$data =~ s/^,// if defined $data;
	return($data);
}

sub bismark_genome_preparation {
	my ($geneIndexFa, $bismark_geneIndexDir, $bowtieOpt, $outLog) = @_;

	my $bismark_genome_preparation_cmd = "bismark_genome_preparation --bowtie2 $bismark_geneIndexDir > $bismark_geneIndexDir/LOG.txt 2>&1";

	LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: Running bismark_genome_preparation\n");
	print_cmd($bismark_genome_preparation_cmd, $outLog);
	#LOG($outLog, "$LCY\n\t$bismark_genome_preparation_cmd\n$N\n");
	my ($geneIndexFaMD5, $temp, $geneIndexFaTempMD5File)  = getMD5($geneIndexFa);
	my $oldMD5File = "$bismark_geneIndexDir/Bisulfite_Genome/MD5SUM";
	#die "geneIndexFaMD5=$LCY$geneIndexFaTempMD5File$N\n";

	my $run_boolean = "\t${LGN}WAS NOT RUN (Using previously made files)$N:";
	my ($check, $md5sum, $md5sum2) = (0,$geneIndexFaMD5,"md5sum2");
	my $bismark_geneIndexDir_exist = (-d "$bismark_geneIndexDir/Bisulfite_Genome/" and -e $geneIndexFaTempMD5File) ? 1 : 0;


	if ($bismark_geneIndexDir_exist == 1 and not defined $force{1}) { # in case in the future we implement -G for bismark folder
		LOG($outLog, "\t#Older bismark folder Bisulfite_Genome $CY$bismark_geneIndexDir/Bisulfite_Genome/$N exist! Checking MD5 if they're the same as current fasta file.\n","NA");
		($md5sum)  = `cat $oldMD5File` =~ /^(\w+)($|[\t ].+)$/;
		($md5sum2) = getMD5($geneIndexFa);
		$bismark_geneIndexDir_exist = 0 if not defined $md5sum or not defined $md5sum2;
		if (defined $md5sum and defined $md5sum2) {
			$bismark_geneIndexDir_exist = 0 if $md5sum ne $md5sum2;
		}
	}
	$run_boolean = "\t${LGN}RAN$N:" if $bismark_geneIndexDir_exist == 0 or defined $force{1};

	#die "bismark folder exist = $bismark_geneIndexDir_exist\n$LCY$bismark_geneIndexDir$N\n$geneIndexFaMD5\n\n";
	if ($bismark_geneIndexDir_exist == 0 or defined $force{1}) {
		LOG($outLog, "\t#Either bismark folder didn't exist or older bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n") if defined $md5sum and $bismark_geneIndexDir_exist == 0;
		LOG($outLog, "\t#Rerunning bismark_genome_preparation due to user request (-F 1)\n") if defined $force{1};
		system($bismark_genome_preparation_cmd) == 0 or die "Failed to run bismark genome preparation: $!\n";
	  # if ($md5script eq "md5") {
	  #    my ($md5res) = `$md5script $geneIndexFa` =~ /^.+\= (\w+)$/; die "Failed to $md5script $geneIndexFa: $!\n" if not defined $md5res;
	  #    system("echo '$md5res $geneIndexFa' > $bismark_geneIndexDir/Bisulfite_Genome/.$md5script") == 0 or LOG($outLog, "Failed to $md5script $geneIndexFa: $!\n");
	  # }
		LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: ${GN}SUCCESS$N: $CY$bismark_geneIndexDir\/Bisulfite_Genome$N made\n");
	}
	else { 
		LOG($outLog, "\n" . date() . "${LPR}bismark_genome_preparation .fa .geneIndex/$N: ${GN}SUCCESS$N: $CY$bismark_geneIndexDir\/Bisulfite_Genome$N already exist (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n");
	}
	system("echo $geneIndexFaMD5 > $bismark_geneIndexDir/Bisulfite_Genome/MD5SUM") == 0 or die "Failed to generate $LCY$bismark_geneIndexDir/Bisulfite_Genome/MD5SUM/$N: $!\n";
	#LOG($outLog, "${run_boolean}$YW ::: $bismark_genome_preparation_cmd :::$N\n");
}

sub run_bismark {
	#($BAMFile) = run_bismark($readFile, $outDir, $BAMFile, $opt_F, $outReadLog, $outLog);
	my ($readFile, $outDir, $BAMFile, $seqFile, $geneIndexFile, $SEQ, $outReadLog, $outLog) = @_;
	my $mapLog = "";
	#LOG($outLog, "Running bismark\n");
	$BAMFile =~ s/(.fq|.fastq|.fq.gz|.fastq.gz)_bismark_bt2/_bismark_bt2/;
	my $ext = "bam";
	my ($BAMFilename) = getFilename($BAMFile, "full");
	
	my $run_boolean = "\n\t${LGN}WAS NOT RUN$N:${YW} ";

	#if (-e $BAMFile and not -e "$outDir/$BAMFilename") {
	#	system("/bin/ln -s $BAMFile $outDir/$BAMFilename") == 0 or LOG($outLog, "Failed to /bin/ln $BAMFile $outDir/$BAMFilename: $!\n") and exit 1;
	#}

	if (not defined $force{2} and -e $BAMFile and not defined $opt_p) {
		#DIELOG($outLog, "bad debug\n") if defined $opt_0;
		my ($BAMFilesize) = `/bin/ls -s $BAMFile` =~ /^(\d+)($|[ \t]+.+)$/;
		$BAMFilesize = 0 if not defined $BAMFilesize;
		my $approx_linecount = int($BAMFilesize / 4+0.5);
		LOG($outLog, "\n" . date() . " ${GN}SUCCESS$N: Output already exist: $CY$BAMFile$N (size=$LGN$BAMFilesize$N, approx_linecount=$LGN$approx_linecount$N lines)\n");
		my $run_bismark_cmd       = "bismark --non_directional -o $outDir --temp_dir $outDir --ambig_bam --ambiguous --unmapped $bowtieOpt $bismark_geneIndexDir $readFile > $outDir/.bismark_log 2>&1";
		$run_boolean .= "${YW}::: $run_bismark_cmd :::$N";
		my ($bismark_report) = "$BAMFile" =~ /^(.+).$ext/; $bismark_report .= "_SE_report.txt";
		if (defined $bismark_report and -e $bismark_report) {
			my ($header, $report, $mapLogTEMP) = parse_bismark_report($bismark_report, $outLog);
			$mapLog = $mapLogTEMP;
			LOG($outLog, "\t#MAP_RESULT\n$LGN$mapLog$N\n");
		}
		else {
			($bismark_report) = $BAMFile =~ /^(.+_bismark_bt2).$ext/; $bismark_report .= "_footLoop_report.txt";
		   ($bismark_report) = <$outDir/*_footLoop_report.txt> if not -e $bismark_report;
		}
		if (not -e $bismark_report) {
			LOG($outLog, "footLoop.pl::run_bismark: can't find bismark_report _footLoop_report.txt $LCY$bismark_report$N!\n");
			$mapLog = "cannot_find_bismark_report";
		}
		else {
			LOG($outLog, "\t#Found parallel-ran bismark report $LCY$bismark_report$N\n");
			my ($header, $report, $mapLogTEMP) = parse_bismark_report($bismark_report, $outLog,1);
			$mapLog = $mapLogTEMP;
		}
		LOG($outLog, "\t#MAP_RESULT\n$LGN$mapLog$N\n");
	}
	else {
		#LOG($outLog, "\t#Rerunning bismark due to user request (-F 2)\n") if defined $force{2};
		$run_boolean = "\n$YW\t ";

		if (not defined $opt_p) { 		# not parallel
			DIELOG($outLog, "BAD DEBUG NOT PARALEL\n");
			DIELOG($outLog, "bad debug\n") if defined $opt_0;
			LOG($outLog, "\t#Not parallel\n");
			my $run_bismark_cmd       = "bismark --non_directional -o $outDir --temp_dir $outDir --ambig_bam --ambiguous --unmapped $bowtieOpt $bismark_geneIndexDir $readFile > $outDir/.bismark_log 2>&1";
			my $run_bismark_cmd_print = "bismark --non_directional -o \$outDir --temp_dir \$outDir --ambig_bam --ambiguous --unmapped \$outDir \$bowtieOpt \$bismark_geneIndexDir \$readFile > \$outDir/.bismark_log 2>&1";
			#my $run_bismark_cmd       = "bismark --non_directional -o $outDir $bowtieOpt $bismark_geneIndexDir $readFile > $outDir/.bismark_log 2>&1";
			#my $run_bismark_cmd_print = "bismark --non_directional -o \$outDir \$bowtieOpt \$bismark_geneIndexDir \$readFile > \$outDir/.bismark_log 2>&1";
			$run_boolean .= " ::: $run_bismark_cmd :::";
			LOG($outLog, "$LCY\n\t$run_bismark_cmd\n$N\n");
			LOG($outReadLog, "footLoop.pl,bismark,$run_bismark_cmd\n","NA");

			my $result = system($run_bismark_cmd);

			if ($result != 0) {
				LOG($outLog, "\t\t${LRD}bismark failed!$N In case bisulfte_genome is corrupted, we're re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_geneIndexDir\n");
				LOG($outReadLog, "$LCY\n\tfootLoop.pl,bismark,bismark_genome_preparation --bowtie2 $bismark_geneIndexDir\n$N","NA");
				system("bismark_genome_preparation --bowtie2 $bismark_geneIndexDir") == 0 or die "Failed to run bismark genome preparation: $!\n";
				LOG($outReadLog, "footLoop.pl,bismark,$run_bismark_cmd");
				system($run_bismark_cmd) == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
			}
			LOG($outLog, "\n" . date() . " ${GN}SUCCESS$N: Output $BAMFile\n");
			my ($bismark_report) = $BAMFile =~ /^(.+_bismark_bt2).$ext/; $bismark_report .= "_SE_report.txt";
			   ($bismark_report) = <$outDir/*bt2_SE_report.txt> if not -e $bismark_report;
			if (not defined $bismark_report or (defined $bismark_report and not -e $bismark_report)) {
				LOG($outLog, "footLoop.pl::run_bismark: can't find bismark_report _SE_report.txt $LCY$bismark_report$N!\n");
			}
			else {
				my ($header, $report, $mapLogTEMP) = parse_bismark_report($bismark_report, $outLog);
				$mapLog = $mapLogTEMP;
			}
			if (not -e $BAMFile) {
				DIELOG($outLog, "footLoop.pl::run_bismark: Fatal error! bam file not found: $LCY$BAMFile$N\n");
			}
			LOG($outLog, "\t#MAP_RESULT\n$LGN$mapLog$N\n");
		}
		else { #parallel

			############
			#  PARALEL #
			############

			my $outFolder = $outDir . "/.run_bismark/";
			system("mkdir -p $outFolder") if not -d $outFolder;
			

			if (defined $force{11} and -d $outFolder) {
				#die "ASDF\n";;
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
			LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: Running ${LCY}footLoop_SplitFastq.pl$N to split readFile into $LGN$parallel$N reads per file\n"); 
			my $splitFastq_resultFile = "$outFolder/.footLoop_SplitFastq_result.txt";
			my $splitFastq_cmd = $parallel ne -1 ? "$footLoop_script_folder/footLoop_SplitFastq.pl -i $readFile -o $outFolder -n $parallel > $splitFastq_resultFile" :
																"/bin/ln -s $readFile $outFolder/$readFilename.0.part.gz > $splitFastq_resultFile";
			print_cmd($splitFastq_cmd, $outLog);

			LOG($outReadLog, "footLoop.pl,run_bismark,$splitFastq_cmd\n","NA");

			my @fastqFiles = <$outFolder/*.part.gz>;
			my $totalfastqFiles = @fastqFiles;
			if (defined $force{3} or not -e $splitFastq_resultFile or $totalfastqFiles == 0) {
				#die "ASDF\n";;
				if (not defined $opt_0) {
					#DIELOG($outLog, "BAD DEBUG\n");
					system($splitFastq_cmd) == 0 or DIELOG($outLog, "footLoop.pl::run_bismark: Failed to run splitFastq_cmd: $LRD$!$N\n$LCY\n$splitFastq_cmd\n$N\n");
				}
				@fastqFiles = <$outFolder/*.part.gz>;
				$totalfastqFiles = @fastqFiles;
				LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: ${GN}SUCCESS$N: Fastq was split successfully into $LGN$totalfastqFiles$N .part.gz files in folder $LCY$outFolder/*.part.gz$N\n");
			}
			else {
				LOG($outLog, "\n" . date() . "${LPR}SplitFastq .fq.gz .part.gz$N: ${GN}SUCCESS$N: Using previously run footLoop_SplitFastq.pl $LGN$totalfastqFiles$N .part.gz files in folder $LCY$outFolder/*.part.gz$N\n");
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

			my @bamFiles = <$outFolder/*_bismark/*bismark_bt2.bam>;#.$ext>; 
			#my @bamFiles      = <$outFolder/*_bismark/*.$ext>;
			my $totalbamFiles = @bamFiles;

			if (defined $force{4} or $totalfastqFiles > $totalbamFiles or $totalbamFiles == 0) {# or not -e $BAMFile) {
				#die "ASDF\n";;
				sbatch_these($bismark_cmd, "bismark", "bam", \@fastqFiles, $max_slurm_job, $outLog);
				#@bamFiles      = <$outFolder/*_bismark/*.$ext>;
				@bamFiles = <$outFolder/*_bismark/*bismark_bt2.bam>;#.$ext>; 
				$totalbamFiles = @bamFiles;
				LOG($outLog, "\n" . date() . "${LPR}bismark .part.gz .bam$N: ${GN}SUCCESS$N: Created $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$outFolder/*_bismark/$N\n");
			}
			else {
				LOG($outLog, "\n" . date() . "${LPR}bismark .part.gz .bam$N: ${GN}SUCCESS$N: Using previously run $LGN$totalbamFiles$N .bam files from bismark run in folders: $LCY$outFolder/*_bismark/$N\n");
			}
			#bismark_check_fq_bam(\@fastqFiles, \@bamFiles, $outLog);
			LOG($outLog, "==================================================\n");			


			################
			# FIX BAM FILE #
			################
			LOG($outLog, "\n" . date() . "${LPR}fixBAM .bam .fixed.gz$N: Running ${LCY}footLoop_2_fixBAMFile.pl$N in parallel\n");

			my $fix_BAMFile_cmd = "footLoop_2_fixBAMFile.pl -b FILENAME -g $seqFile -i $geneIndexFile -o FOLDER";
			print_cmd("$fix_BAMFile_cmd", $outLog);

			my @fixedFiles = <$outFolder/*_bismark/*.$ext.fixed.gz>;
			my $totalfixedFiles = @fixedFiles;

			if (defined $force{5} or $totalfixedFiles < $totalbamFiles) {
				#die "ASDF\n";;
				sbatch_these($fix_BAMFile_cmd, "fixBAM", "fixed.gz", \@bamFiles, $max_slurm_job, $outLog);
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
			LOG($outLog, "\n" . date() . "${LPR}filterBAM .fixed.gz .filtered.gz$N: Running ${LCY}footLoop_3_filterBAMFile.pl$N in parallel\n");

			my $filter_BAMFile_cmd = "footLoop_3_filterBAMFile.pl -r $opt_r -b FILENAME -g $seqFile -L $opt_L -O FOLDER -o FOLDER -i $geneIndexFile2";
			   $filter_BAMFile_cmd .= " -q $opt_q" if defined $opt_q;
			   $filter_BAMFile_cmd .= " -x $opt_x" if defined $opt_x;
			   $filter_BAMFile_cmd .= " -y $opt_y" if defined $opt_y;
			print_cmd("$filter_BAMFile_cmd", $outLog);

			my @filterBAMFile = <$outFolder/*_bismark/*.log.txt>; #log.txt coz there could be massive number of .filtered files which will slow down the script
			my $totalfilterBAMFile = @filterBAMFile;

			if (defined $force{6} or @fixedFiles < @bamFiles or $totalfilterBAMFile == 0) {
				sbatch_these($filter_BAMFile_cmd, "filterBAM", "log.txt", \@bamFiles, $max_slurm_job, $outLog);
				@filterBAMFile = <$outFolder/*_bismark/*.log.txt>; #\_filterBAM/*.log.txt>;
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
					#die "ASDF\n";;
					#DIELOG($outLog, "BAD DEBUG\n");
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
							#DIELOG($outLog, "BAD DEBUG\n");
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
					#die "ASDF\n";;
					#DIELOG($outLog, "BAD DEBUG\n");
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
					LOG($outLog, " $q: $LCY$header[$q]$N=$LGN$report[$q]$N","NA") if $p < 5;
					LOG($outLog, "\nfootLoop.pl::run_bismark: parse_bismark_report loop, ERROR undefined header at i=$q\n") if not defined $header[$q];
					LOG($outLog, "\nfootLoop.pl::run_bismark: parse_bismark_report loop, ERROR header ($header[$q]) isnt' same as HEADER ($HEADER[$q]) at i=$q\n") if defined $header[$q] and $header[$q] ne $HEADER[$q];
					# just next coz we're logging, not that important
					next if not defined $header[$q];
					next if $header[$q] ne $HEADER[$q];
					next if not defined $report[0];
					next if not defined $report[$q];
					$REPORT[$q] += $report[$q] if $header[$q] !~ /^perc_/;
					$REPORT[$q] += $report[0]*$report[$q]/100 if $header[$q] =~ /^perc_/; #1000*98.5/100 = 985 * 2600 = 2,500,000
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
				$REPORT[$q] = int(10000*$REPORT[$q]/($REPORT[0])+0.5)/100 if $HEADER[$q] =~ /^perc_/; # 2,500,000/(1000*2600) = 98.5
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
			LOG($outLog, "==================================================\n");			

		}
	}
	#LOG($outLog, "${run_boolean}$N\n");#::: bismark $bowtieOpt $bismark_geneIndexDir $readFile :::$N\n");

	# print to $outReadLog
	LOG($outReadLog, "footLoop.pl,BAMFile,$outDir/$BAMFilename\n","NA");
	LOG($outReadLog, $mapLog,"NA");

	return("$outDir/$BAMFilename");
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


#sub run_bismark_methylation_extractor {
#	my ($BAMFileInput, $max_slurm_job, $outLog) = @_;
#	my @BAMFileArr = $BAMFileInput =~ /ARRAY/ ? @{$BAMFileInput} : ($BAMFileInput);
#	#my $outDir = "$BAMFile\_bismark_methylation_extractor";
#	#system("mkdir -p $outDir") == 0 or DIELOG($outLog, "\n" . date() . " Failed to mkdir -p $LCY$outDir$N: $!\n\n");
#	my $bme_cmd = "
#		bismark_methylation_extractor --version --gzip --comprehensive --merge_non_CpG \
#		--output FILENAME\_bismark_methylation_extractor \
#		--split_by_chromosome
#	";
#	LOG($outLog, "\n" . date() . "2c. Running bismark_methylation_extractor\n\n");
#	LOG($outLog, "$LCY$bme_cmd$N\n\n");
#	sbatch_these($bme_cmd, \@BAMFileArr, $max_slurm_job, $outLog);
#}

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

sub sbatch_these {
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

		#if ($suffix eq "fixBAM") {
		#	#LOG($outLog, "$file.fixed.gz\n");
		#	if (-e "$file\.fixed.gz") {#not -e "$file\_filterBAM.done" and not -e "$file\_fixBAM.done" and -e "$file\_logBAM.done") {
		#		system("echo DONE > $file\_fixBAM.done");
		#	#	if (-e "$file\_logBAM.done" and -e "$file\_logBAM.sbatch" and -e "$file\_logBAM.sbout") {
		#			system("/bin/mv $file\_filterBAM.sbatch $file\_fixBAM.sbatch") if -e "$file\_filterBAM.sbatch" and not -e "$file\_fixBAM.sbatch";
		#			system("/bin/mv $file\_filterBAM.sbout $file\_fixBAM.sbout") if -e "$file\_filterBAM.sbout" and not -e "$file\_fixBAM.sbout";
		#			system("/bin/mv $file\_filterBAM.done $file\_fixBAM.done") if -e "$file\_filterBAM.done" and not -e "$file\_fixBAM.done";
		#	#	}
		#		#DIELOG($outLog, "GOOD DEBUG\n") if $i == 0;
		#	}
		#}
			#if (-e "$file\_logBAM.done" and -e "$file\_logBAM.sbatch" and -e "$file\_logBAM.sbout") {
			#	LOG($outLog, "${LGN}DOING$N /bin/mv $file\_filterBAM.sbatch $file\_fixBAM.sbatch\n");
			#	#DIELOG($outLog, "GOOD EBUG\n");
			#	system("/bin/mv $file\_filterBAM.sbatch $file\_fixBAM.sbatch") if -e "$file\_filterBAM.sbatch";
			#	system("/bin/mv $file\_filterBAM.sbout $file\_fixBAM.sbout") if -e "$file\_filterBAM.sbout";
			#	system("/bin/mv $file\_filterBAM.done $file\_fixBAM.done") if -e "$file\_filterBAM.done";
			#	exit 0;
			#}
			#else {
			#	LOG($outLog, "${LRD}NOT DOING$N /bin/mv $file\_filterBAM.sbatch $file\_fixBAM.sbatch\n");
			#}
			#DIELOG($outLog, "GOOD EBUG\n") if $i == 9;
		#}
		#if ($suffix eq "filterBAM" and -e "$file\_logBAM.done") {
		#	#DIELOG($outLog, "GOOD EBUG\n");
		#	if (-e "$file\.fixed.gz") {
		#		#LOG($outLog, "${LGN}DOING$N /bin/mv $file\_logBAM.sbatch $file\_filterBAM.sbatch\n");
		#		system("/bin/mv $file\_logBAM.sbatch $file\_filterBAM.sbatch") if -e "$file\_logBAM.sbatch";
		#		system("/bin/mv $file\_logBAM.sbout $file\_filterBAM.sbout") if -e "$file\_logBAM.sbout";
		#		system("/bin/mv $file\_logBAM.done $file\_filterBAM.done") if -e "$file\_logBAM.done";
		#	}
		#	#else {
		#	#	#LOG($outLog, "${LRD}NOT DOING$N /bin/mv $file\_logBAM.sbatch $file\_filterBAM.sbatch\n");
		#	#}
		#	#DIELOG($outLog, "GOod debug\n");
		#}
		
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
			$sbatchprint .= "conda activate footLoop2\n";
			$sbatchprint .= "$cmdcopy && echo \"Done!\" > $donefile\n\n";

		open (my $out, ">", $sbatchfile) or die "Can't write to $LCY$sbatchfile$N: $!\n";
		print $out $sbatchprint;
		close $out;

		if (not defined $force{0} and not defined $force_sbatch and -e $donefile) {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles sbatch_these $suffix$N: sbatch $LCY$sbatchfile$N # ${LGN}DONE$N\n");
			next;
		}
		else {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles sbatch_these $suffix$N: sbatch: $LCY$sbatchfile$N\n");
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
		#system("touch $file.done") == 0 or die "failed to touch $file.done: $!\n";
	}
	my $sleep = 0;
	#while (1) {
	#	my ($job_left) = squeue_check($jobidhash);
	#	LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 60 == 0;
	#	last if ($job_left < $max_parallel_run);
	#	$sleep ++;
	#	sleep 1;
	#}
	#$sleep = 0;
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
	#my $mapLog  = "footLoop.pl,map," . "header\tlabel\t" . join("\t", @header) . "\tfootLoop_outDir\tuuid\n";
	#   $mapLog .= "footLoop.pl,map," . "record\t$label\t" . join("\t", @report) . "\t$outDir\t$uuid\n";
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
		#LOG($outLog, "$YW\t::: /bin/cp $geneIndexFile $geneIndexFileNew :::$N\n");
		my $buffer_cmd = "/bin/cp $geneIndexFile $geneIndexFileNew";
		print_cmd($buffer_cmd, $outLog);
		system($buffer_cmd) == 0 or DIELOG($outLog, "\tbedtools_bed_change: Failed because of $LGN$!$N cmd=$LCY$buffer_cmd$N\n\n");
		LOG($outLog, "\n" . date() . "${LPR}bedtools_bed_change .bed _x_y.bed$N: ${GN}SUCCESS$N: Copied ${LCY}geneIndexFile$N $LCY$geneIndexFileNew$N\n");
	}
	else {
		my $buffer_cmd = "$footLoop_script_folder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1";
		print_cmd($buffer_cmd, $outLog);
		system($buffer_cmd) == 0 or DIELOG($outLog, "\tbedtools_bed_change: Failed because of $LGN$!$N cmd=$LCY$buffer_cmd$N\n\n");
		#system($buffer_cmd) == 0 or DIELOG($outLog, "Failed to get (beg $bufferL, end $bufferR) bp of $geneIndexFile!\n");
		LOG($outLog, "\n" . date() . "${LPR}bedtools_bed_change _x_y.bed$N: ${GN}SUCCESS$N: Created ${LCY}geneInexFile$N $LCY$geneIndexFileNew$N\n");
	}
	LOG($outLog, "==================================================\n");			

	return($geneIndexFileNew);
}
#system("$footLoop_script_folder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1") == 0 or LOG($outLog, "\tfootLoop.pl::get_geneIndex_fasta: Failed to$YW $footLoop_script_folder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1$N\n: $!\n") and exit 1;

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

	#LOG($outLog, "\n\tb. Getting fasta sequence from bed file$LCY $geneIndexFile$N genome file$LCY $genomeFile$N\n");

	LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: Getting fasta sequence using ${LCY}fastaFromBed$N\n");
	my $fastaFromBed_cmd = "fastaFromBed -fi $genomeFile -bed $geneIndexFile -fo $geneIndexFaTemp -nameOnly";
	print_cmd($fastaFromBed_cmd, $outLog);
	#LOG($outLog, "$YW\t::: $cmd :::$N\n");
	my $fastaFromBed_cmd_success = system("$fastaFromBed_cmd");
	if ($fastaFromBed_cmd_success == 0) {
		LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: ${GN}SUCCESS$N: Created fasta file $LCY$geneIndexFaTemp$N\n");
	}
	else {
		LOG($outLog, "\n" . date() . "${LPR}fastaFromBed .bed .fa$N: ${LRD}FAIL$N: Failed to create fasta file $LCY$geneIndexFaTemp$N: $LGN$!$N\n");
	}
	LOG($outLog, "==================================================\n");			

	#$geneIndexFaTemp = uppercaseFasta($geneIndexFaTemp, $outLog);
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
	#my ($SEQ, $geneIndexHash, $seqFile, $bismark_geneIndexDir) = parse_geneIndexFile($geneIndexFile, $outDir, $outLog, $seqFile, $minReadL);
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

	# -r read not defined
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
			#my $prev_opt_g = $opt_g;
	
			my $curr_opt_g = "$opt_n/$faFile";
			my $faiFile    = "$opt_n/$faFile.fai";
			$opt_i         = "$opt_n/$faFile.fai.bed";
	
			if (not -e $curr_opt_g) {#"$opt_n/$faFile") {
				my ($ln_success) = system("/bin/ln -s $prev_opt_g $opt_n/$faFile");# == 0 or die "footLoop.pl::check_sanity: can't create symlink of fasta file into output folder: $!\n$LCY$opt_g$N\n$LGN$opt_n$N\n\n";
				#print "ln_success = $LGN$ln_success$N\n";
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
5: filter bam (footLoop_2_fixBAMFile.pl)
6: log bam (footLoop_3_filterBAMFile.pl)
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
./footLoop.pl -r PCB190425.fq.gz -n PCB190425_MAP -g hg19.fa -l PCB190425 -x -10 -y 10 -i geneIndexes.bed
./footPeak.pl -n PCB190425_MAP -o PCB190425_PEAK
./footClust.pl -n PCB190425_PEAK
./footPeak_graph.pl -n PCB190425_PEAK
./footPeakGTF.pl -n PCB190425_PEAK
./footPeak_GCprofile.pl -n PCB190425_PEAK -i geneIndexes.bed
./footStats.pl -n PCB190425_PEAK

";

	return($usageshort, $usage, $usagelong, $example);
}


__END__
my $bismarkloc  = `which bismark` ; chomp($bismarkloc) ; $bismarkloc  = "N/A" if $bismarkloc eq "";
my $bedtoolsloc = `which bedtools`; chomp($bedtoolsloc); $bedtoolsloc = "N/A" if $bedtoolsloc eq "";
my $bowtie2loc  = `which bowtie2` ; chomp($bowtie2loc) ; $bowtie2loc  = "N/A" if $bowtie2loc eq "";
my $samtoolsloc = `which samtools`; chomp($samtoolsloc); $samtoolsloc = "N/A" if $samtoolsloc eq "";
my $Rscriptloc  = `which Rscript` ; chomp($Rscriptloc) ; $Rscriptloc  = "N/A" if $Rscriptloc eq "";
my $bismarktest  = `bismark  --version 2>&1 | head`; chomp($bismarktest);
my $bedtoolstest = `bedtools --version 2>&1 | head`; chomp($bedtoolstest);
my $bowtie2test  = `bowtie2  --version 2>&1 | head`; chomp($bowtie2test);
my $samtoolstest = `samtools --version 2>&1 | head`; chomp($samtoolstest);
my $Rscripttest  = `Rscript  --version 2>&1 | head`; chomp($Rscripttest);
my $print = "

# Locations:
- bedtools = $LCY$bedtoolsloc$N
- samtools = $LCY$samtoolsloc$N
- bowtie2 = $LCY$bowtie2loc$N
- bismark = $LCY$bismarkloc$N
- R = $LCY$Rscriptloc$N

# Versions:
bedtools --version:
${LGN}>>$bedtoolstest<<$N

samtools --version:
${LGN}>>$samtoolstest<<$N

bowtie2 --version:
${LGN}>>$bowtie2test<<$N

bismark --version:
${LGN}>>$bismarktest<<$N

Rscript --version:
${LGN}>>$Rscripttest<<$N

print $print . "\n";
print "\n" if $print =~ /not found/i;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoop_script_folder = dirname(dirname abs_path $0);
$footLoop_script_folder .= "/footLoop";
my $cmd = "$footLoop_script_folder/check_software.pl | tail -n 12";
my @version = `$cmd`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}


__END__
LABEL
elsif ($outDir =~ /PCB\d+/i) {
   ($label) = $outDir =~ /(PCB\d+)/i;
   $label = uc($label);
   $label =~ s/PCB0(\d+)/PCB$1/ if $label =~ /PCB0\d+/;
}
if (not defined $label and $outDir =~ /PCB[0\-_]*\d+/i) {
   ($label) = $outDir =~ /(PCB[0\-_]*\d+)/i;
   $label = uc($label);
   $label =~ s/PCB[0\-_]+(\d+)/PCB$1/g if $label =~ /PCB[0\-_]+\d+/;
}
if (not defined $label and $label =~ /PCB/) {
   ($label) = $outDir =~ /(PCB.{1,5})\/?/;
}
if (not defined $label) {
   die "Please make sure your output folder (-n) contain PCB(number) e.g. PCB12: -n 180202_PCB12_footLoop_output\n\n";
}

if (not defined $opt_l and $outDir =~ /_BC\d+_PFC\d+_\w+\./i) {
   my ($label2) = $outDir =~ /_(BC\w+)\./i;
   $label = "$label\_$label2";
   $label = uc($label);
   $label =~ s/PCB[0\-_]+(\d+)/PCB$1/g;
}
	$OPTS{l} = $label;

if ($outDir =~ /PCB\d+_bcBC\d+_.+desc\w+/) {
	($label) = $outDir =~ /(PCB\d+_bcBC\d+_.+desc[a-zA-Z0-9]+)_?/;
}
print date() . "FATAL ERROR ON LABEL -l\n" and die if not defined $label;
=cut

__END__
LABEL sanity_chck
=comment
	$errors .= "-n $LGN<output directory>$N is not defined.\n" if not defined($opt_n);
	my $message = "";
	if (defined $opt_n and not -d $opt_n) {
		$message = `mkdir $opt_n 2>&1`;
		chomp($message);
	}
	$errors .= "-n $LGN<output directory>$N is defined ($YW$opt_n$N) but folder does not exist and cannot be created.\n   -> mkdir $YW$opt_n$N returned $LCY$message$N.\n" if defined $opt_n and not -d $opt_n;

	# label error
	my $label;
	if (defined $opt_l) {
		$label = $opt_l;
	}
	elsif (defined $opt_n) {
		my ($readName) = getFilename($opt_n, "full");
		if ($readName =~ /PCB\d+_bcBC\d+_plasmid.+_+desc.+/i) {
			my ($PCB, $BC, $plasmid, $desc) = $readName =~ /PCB(\d+)_bcBC(\d+)_plasmid(.+)_desc([a-zA-Z0-9_]+)\.?/;
			if (defined $desc) {
				$label = "PCB$PCB\_bcBC$BC\_plasmid$plasmid\_desc$desc";
			}
		}
		elsif ($readName =~ /PCB\d+/i) {
			my ($PCB) = $opt_n =~ /PCB(\d+)/;
			if (defined $PCB) {
				$label = "PCB$PCB";
			}
		}
	}
	if (((defined $label and $label !~ /^PCB\d+(_bcBC\d+_plasmid.+_desc.+)$/) or not defined $label) and defined $opt_r) {
		my ($readName) = getFilename($opt_r, "full");
		if ($readName =~ /PCB\d+_bcBC\d+_plasmid.+_+desc.+/i) {
			my ($PCB, $BC, $plasmid, $desc) = $readName =~ /PCB(\d+)_bcBC(\d+)_plasmid(.+)_desc([a-zA-Z0-9_]+)\.?/;
			if (defined $desc) {
				$label = "PCB$PCB\_bcBC$BC\_plasmid$plasmid\_desc$desc";
			}
		}
		elsif ($readName =~ /PCB\d+/i) {
			my ($PCB) = $opt_r =~ /PCB(\d+)/;
			if (defined $PCB) {
				$label = "PCB$PCB";
			}
		}
	}
	$opt_l = $label if not defined $opt_l and defined $label;
	$errors .= "-l $LGN<label>$N: is not defined!\n\t-> $YW Has to be this exact format$N:.\n\t  - Invivo: use PCB<digits> e.g. PCB1 or PCB12636.\n\t  - Invitro: Use PCB<digits>_bcBC<digits>_plasmid<string>_desc<string> e.g. PCB13_bcBC13_plasmidPFC53_descSUPERCOILEDCIRC.\n" if not defined($opt_l);

	# label exists but wrong format
	$errors .= "-l $LGN<label$N: is defined ($YW$opt_l$N) but wrong format (${LPR}case sensitive$N)! PCB<digits> or PCB<digits>_bcBC<digits>_plasmid<string>_desc<string>.\n" if defined $opt_l and $opt_l !~ /^PCB\d+$/ and $opt_l !~ /^PCB\d+_bcBC\d+_plasmid[0-9a-zA-Z]+_desc[0-9a-zA-Z]+$/;
=cut



__END__
			#my $result = system("run_script_in_parallel2.pl -v \"srun -p high --mem 8000 bismark -o $outDir/.run_bismark/ $bowtieOpt $bismark_geneIndexDir FILENAME >> FILENAME.bismark.log 2>&1\" $outFolder .part 20");
			#LOG($outReadLog, "footLoop.pl,run_bismark,\"run_script_in_parallel2.pl -v \\\"srun -p high --mem 8000 bismark -o $outDir/.run_bismark/ $bowtieOpt $bismark_geneIndexDir FILENAME >> FILENAME.bismark.log 2>&1\\\" $outFolder .part 20");




sub print_R_heatmap {
	my ($SEQ) = @_;
	

	my $box;
	foreach my $gene (sort keys %{$SEQ}) {
		my $Rscript = "$filteredDir/.$gene\_MakeHeatmap.R";
		open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
	
		# Breaks and color for R script
		# 0 is bad or not data (6)
		# 1 is non C/G
		# 4 is CH non conv
		# 5 is CG non conv
		# 6 is CH conv
		# 7 is CG conv
		# 8 is CH peak
		# 9 is CG peak
		
		my $breaks = "c(-1"	; my $colors = "c("               ; #BEGIN
		   $breaks .= ",0"	;    $colors .= "\"grey\""        ; # 0 is bad or no data (used to be 6)
		   $breaks .= ",1"	;    $colors .= ",\"white\""		 ; # 1 is non CG
		   $breaks .= ",4"	;    $colors .= ",\"cornsilk\""	 ; # 4 is CH non conv
			$breaks .= ",5"	;    $colors .= ",\"peachpuff\""  ; # 5 is CG non conv
			$breaks .= ",6"	;    $colors .= ",\"green4\""  	 ; # 6 is CH conv
			$breaks .= ",7"	;    $colors .= ",\"seagreen4\""  ; # 7 is CG conv
			$breaks .= ",8"	;    $colors .= ",\"red4\""    	 ; # 8 = PEAK Converted CH
			$breaks .= ",9"	;    $colors .= ",\"maroon4\""    ; # 9 = PEAK Converted CG
		   $breaks .= ")"    ;    $colors .= ")"               ; #END

		#Function to turn @seq into number and array for R
		my $seqR = "seq = c(";
		my @seq = @{$SEQ->{$gene}{seq}};
		for (my $s = 0; $s < @seq; $s++) {
			my $nuc = $seq[$s];
			#$nuc = $nuc =~ /A/i ? 11 : $nuc =~ /C/i ? 12 : $nuc =~ /G/i ? 13 : $nuc =~ /T/i ? 14 : 15;
			$nuc = $nuc =~ /A/i ? "seagreen4" : $nuc =~ /C/i ? "blue4" : $nuc =~ /G/i ? "saddlebrown" : $nuc =~ /T/i ? "red4" : "grey";
	
			$seqR .= "\"$nuc\"";
			$seqR .= $s == @seq - 1 ? ")" : ",";
		}
	
		my $labelz .= "p = p + geom_text(aes(x=min(dm\$x),label=id),hjust=0)";# if defined $opt_T;
		my $RBox = ""; my $RBoxPlot = "";
		if (defined $box) {
			if (defined $box->{$gene}) {
				my %RBox;
				foreach my $lines (sort @{$box->{$gene}}) {
					my ($chr, $beg, $end, $geneBox) = split("\t", $lines);
					push(@{$RBox{name}}, $geneBox);
					push(@{$RBox{beg}}, $beg);
					push(@{$RBox{end}}, $end);
					print "\t$geneBox, beg=$beg, end=$end\n";
				}
				if (defined $RBox{name}) {
					my $RBoxname = "RBoxname = c(\"" . join("\",\"", @{$RBox{name}}) . "\")";
					my $RBoxbeg = "RBoxbeg = c(" . join(",", @{$RBox{beg}}) . ")";
					my $RBoxend = "RBoxend = c(" . join(",", @{$RBox{end}}) . ")";
					$RBox = "
						mydm_xmin = min(dm\$x)
						mydm_xmax = max(dm\$x)
						$RBoxname
						$RBoxbeg
						$RBoxend
						if (length(RBoxbeg[RBoxbeg < mydm_xmin]) > 0) {
							RBoxbeg[RBoxbeg < mydm_xmin] = mydm_xmin
						}
						if (length(RBoxbeg[RBoxbeg > mydm_xmax]) > 0) {
							RBoxbeg[RBoxbeg > mydm_xmax] = mydm_xmax
						}
						if (length(RBoxend[RBoxend < mydm_xmin]) > 0) {
							RBoxend[RBoxend < mydm_xmin] = mydm_xmin
						}
						if (length(RBoxend[RBoxend > mydm_xmax]) > 0) {
							RBoxend[RBoxend > mydm_xmax] = mydm_xmax
						}
						RBox = data.frame(name=RBoxname,my_xmin=RBoxbeg,my_xmax=RBoxend, my_ymin=min(dm\$y),my_ymax = max(dm\$y))
					";
			
					$RBoxPlot = "+ geom_rect(data=RBox,aes(x=RBox\$my_xmin,y=RBox\$my_ymin,xmin=RBox\$my_xmin,xmax=RBox\$my_xmax,ymin=RBox\$my_ymin+0.5,ymax=RBox\$my_ymax-0.5),size=2,fill=NA,color=\"black\")";
				}
			}
		}
		print $out "
	
			.libPaths()
			library(\"GMD\")
			library(ggplot2)
			library(Cairo)
			#Sequence
			$seqR
	
	
			# print peak
			files = c()
			PDFout = c()
			PNGout = c()
			PEAKS = paste(files,\".PEAKS\",sep=\"\")
	
			info = file.info(files)
			files = rownames(info[info\$size != 0, ])
			for (i in 1:length(files)) 
			{
				if (file.exists(files[i])) 
				{
					myclust = TRUE
					if(length(grep(\"_NOPK\", files[i])) != 0) {
						myclust = FALSE
					}
					
					df = read.table(files[i],row.names=1,sep=\"\\t\")
					if (i > 2 & dim(df)[1] > 1000) {
						df = df[seq(1,1000),]
					}
	
					# in case there's only 1 peak:
					if (dim(df)[1] == 1) {
						df[2,] = df[1,]
					}
					if (length(df) > 0) {
	
	               rownames(df) = gsub(\"^SEQ_(\\\\w+)\$\",\"\\\\1\",rownames(df),perl=T)
	               ### Main ###
	               library(ggplot2)
	               library(reshape2)
	               library(grid)
	               h = hclust(d=dist(df))
	               df = df[h\$order,]
	               colnames(df) = seq(1,dim(df)[2])
	               df\$y = seq(1,dim(df)[1])
	               df\$id = rownames(df)
	               dm = melt(df,id.vars=c(\"id\",\"y\"))
	               colnames(dm) = c(\"id\",\"y\",\"x\",\"value\")
	               dm\$x = as.numeric(as.character(dm\$x))
	               dm\$y = as.numeric(as.character(dm\$y))
	               dm\$value = as.factor(as.character(dm\$value))
	
						if (file.exists(PEAKS[i])) {
	  	            	mypeaks = read.table(PEAKS[i]);
	   	            colnames(mypeaks) = c(\"id\",\"x\",\"xend\")
	      	         mydf = data.frame(id=rownames(df), y=df\$y)
	         	      mypeaks = merge(mypeaks, mydf, by=\"id\")
						} else {
							mypeaks=NULL
						}
					
	####### R BOX
	
					$RBox
	
	#############
	
	               p = ggplot(dm,aes(x,y)) +
	                   geom_tile(aes(fill=value)) + coord_fixed(ratio=4) +
	                   scale_fill_manual(values=c(\"0\"=\"cornsilk\",\"1\"=\"green4\",\"2\"=\"white\",\"3\"=\"peachpuff\",\"4\"=\"seagreen4\",\"5\"=\"maroon4\",\"6\"=\"grey\",\"9\"=\"red4\")) +
	                   scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	                   theme(line = element_blank(),
	                         axis.text = element_blank(),
	                         axis.title = element_blank()
	                   ) + ggtitle(paste(\"(\",dim(df)[1],\" peaks)\",sep=\"\")) $RBoxPlot
	################# R BOX PLOT
	
						if (length(mypeaks) > 0) {
							p = p + geom_rect(data=mypeaks,aes(xmin=mypeaks\$x,xmax=mypeaks\$xend,ymin=mypeaks\$y-0.5,ymax=mypeaks\$y+0.5),color=rgb(1,0,0),size=0.2,fill=NA); #+
									#geom_text(data=mypeaks,aes(x=mypeaks\$xend, y=mypeaks\$y, label=mypeaks\$id))
						}
						
						if (length(df) < 100) {
							$labelz
						}
	
	               gt <- ggplot_gtable(ggplot_build(p))
	               ge <- subset(gt\$layout, name == \"panel\")
	
	               ### PNG ###
	
	               print(paste(\"    \",i,\"B. Doing PNG of \",files[i],sep=\"\"))
	               png(type=\"cairo\",PNGout[i],height=4*2*dim(df)[1],width=2*dim(df)[2])
	               grid.draw(gt[ge\$t:ge\$b, ge\$l:ge\$r])
	               dev.off()
	
	
	
	               ### PNG ###
	               print(paste(i, \". Doing pdf of \",files[i],sep=\"\"))
	               pdf(PDFout[i])#,height=4*2*dim(df)[1],width=2*dim(df)[2]/4)
	               grid.draw(gt[ge\$t:ge\$b, ge\$l:ge\$r])
	               dev.off()
	               
					}
				}
			}
		";
		close $out;
		my $cmd = "R --no-save < $Rscript >> $logFile\n";
		print ($cmd);
		system($cmd) if -e $Rscript and -s $Rscript > 10;
		system("rm Rplots.pdf") if -e "Rplots.pdf";
		LOG($outLog, "\n${LPR}If there is not PDF made from step (8) then it's due to too low number of read/peak$N\n\n");
	}
}



__END__
#split_BAMFile {
	my ($BAMFile, $seqFile, $outReadLog, $outLog) = @_;
	LOG($outLog, "\n\ta. Fixing BAM file $CY$BAMFile$N with $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl\n");
	my ($BAMMD5) = getMD5($BAMFile);
	my $filteredDir = "$outDir/.BAMFile_$BAMMD5/";
	check_if_result_exist(["$filteredDir/.GOOD"], $outLog);
	makedir("$filteredDir") if not -d "$filteredDir";
	my $checkBAM = 1;
	my ($BAMFileName) = getFilename($BAMFile, "full");
	my $BAMFileGZ = "$filteredDir/$BAMFileName.fixed.gz";
	$checkBAM = 0 if not -e "$filteredDir/$BAMFileName.fixed" and not -e "$filteredDir/$BAMFileName.fixed.gz";
	makedir($filteredDir);
	if (defined $opt_F) {
		$checkBAM = 0;
		if (-e "$filteredDir/$BAMFileName.fixed.gz") {
			my ($backup_log) = create_backup("$filteredDir/$BAMFileName.fixed.gz", "mv");
			LOG($outLog, $backup_log);
		}
		if (-e "$filteredDir/$BAMFileName.fixed" and $checkBAM == 0) {
			my ($backup_log) = create_backup("$filteredDir/$BAMFileName.fixed", "mv");
			LOG($outLog, $backup_log);
		}
	}
	else {
		if (-e "$filteredDir/$BAMFileName.fixed.gz") {
			my ($BAMLineCount2) = linecount("$filteredDir/$BAMFileName.fixed.gz");
			my ($BAMLineCount1) = linecount($BAMFile);
			$checkBAM = $BAMLineCount1 - 10 > $BAMLineCount2 ? 0 : 2;
			if ($checkBAM eq 0) {
				my ($backup_log) = create_backup("$filteredDir/$BAMFileName.fixed.gz", "mv");
				LOG($outLog, $backup_log);
				LOG($outLog, "\tfootLoop.pl subroutine split_BAMFile:: fixed BAM file $LCY$filteredDir/$BAMFileName.fixed.gz$N exists but total row is less than total BAMFile $BAMFile row ($BAMLineCount1 - 500 > BAMFile.fixed.gz: $BAMLineCount2)!\n");
			}
			else {
				LOG($outLog, "\n" . date() . "\tfootLoop.pl subroutine split_BAMFile::$LGN SUCCESS!!$N fixed BAM file $LCY$filteredDir/$BAMFileName.fixed.gz$N exists (MD5=$LGN$BAMMD5$N) and total row $LGN($BAMLineCount2)$N >= total BAMFile row $LGN($BAMLineCount1 - 500)$N ($LCY$BAMFile$N)!\n") if $checkBAM == 2;
			}
		}
		if (-e "$filteredDir/$BAMFileName.fixed" and $checkBAM == 0) {
			my ($BAMLineCount2) = linecount("$filteredDir/$BAMFileName.fixed");
			my ($BAMLineCount1) = linecount($BAMFile);
			$checkBAM = $BAMLineCount1 - 10 > $BAMLineCount2 ? 0 : 1;
			if ($checkBAM eq 0) {
				my ($backup_log) = create_backup("$filteredDir/$BAMFileName.fixed", "mv");
				LOG($outLog, $backup_log);
				LOG($outLog, "\tfootLoop.pl subroutine split_BAMFile:: .gz does not exist and fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N exists but total row is less than total BAMFile $BAMFile row ($BAMLineCount1 - 500 > BAMFile.fixed: $BAMLineCount2)!\n");
			}
			else {
				LOG($outLog, "\n" . date() . "\tfootLoop.pl subroutine split_BAMFile::$LGN SUCCESS!!$N fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N exists (MD5=$LGN$BAMMD5$N) and total row $LGN($BAMLineCount2)$N >= total BAMFile row $LGN($BAMLineCount1 - 500)$N ($LCY$BAMFile$N)!\n");
			}
		}
	}

	if ($checkBAM == 0) {
		LOG($outLog, "\tfootLoop.pl subroutine split_BAMFile:: fixed BAM file $LCY$filteredDir/$BAMFileName.fixed$N or .gz does not exist!\n");
		LOG($outLog, "\t${YW}$footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -s $seqFile -o $filteredDir$N\n");
		system("$footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir") == 0 or LOG($outLog, "Failed to run $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir: $!\n") and exit 1;
		LOG($outReadLog, "footLoop.pl,split_BAMFile,$footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -o $filteredDir\n","NA");
		if (not -e "$filteredDir/$BAMFileName.fixed.gz") {
			LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed");
			system("gzip $filteredDir/$BAMFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $filteredDir/$BAMFileName.fixed: $!\n");
		}
		else {
			LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed\n\t${LGN}Already exist! So not overwriting!\n");
		}
		$checkBAM = 1;
	}
	else {
		LOG($outLog, "\t${LGN}WAS NOT RUN$N: ${YW}::: $footLoop_script_folder/lib/footLoop_2_fixBAMFile.pl -n $outDir -s $seqFile -o $filteredDir :::$N\n");
		LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed\n\t${LGN}Already exist! So not overwriting!\n");
		# rm old (bad) .gz if it exists
		#if (-e $BAMFileGZ) {
		#	my $indice = 1;
		#	my $currbamfilegz = $BAMFileGZ;
		#	while (-e $currbamfilegz) {
		#		$indice ++;
		#		my $currbamfilegz = $BAMFileGZ . ".$indice.gz";
		#	}
		#	my $cmd2 = "mv $BAMFileGZ $currbamfilegz";
		#	system($cmd2) == 0 or die "\nFailed to $LCY$cmd2$N: $!\n\n";
		#	LOG($outLog, "\t{YW}Moving curr fixed.gz to backup!\n\t$LCY$cmd2$N\n");
		#	#LOG($outLog, "\t/bin/rm $BAMFileGZ") if -e '$BAMFileGZ';
		#}
		#system("/bin/rm $BAMFileGZ") == 0 or LOG($outLog, 'Failed to rm $BAMFileGZ: $!\n') if -e '$BAMFileGZ';
		# gzip the new .fixed
		#LOG($outLog, "\tgzip $filteredDir/$BAMFileName.fixed");
		#system("gzip -f $filteredDir/$BAMFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $filteredDir/$BAMFileName.fixed: $!\n");
	}

	# re-md5 BAMfile gz
	LOG($outLog, "\t${YW}$md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5$N\n");
	if ($md5script eq "md5sum") {
		system("$md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	if ($md5script eq "md5") {
		my ($md5res) = `$md5script $BAMFileGZ` =~ /^.+\= (\w+)$/; die "Failed to $md5script $BAMFileGZ: $!\n" if not defined $md5res;
		system("echo '$md5res $BAMFileGZ' > $filteredDir/.$BAMFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $BAMFileGZ > $filteredDir/.$BAMFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	LOG($outReadLog, "footLoop.pl,BAMFixedFile,$BAMFileGZ\n","NA");
	LOG($outReadLog, "footLoop.pl,BAMFixedFileMD5,$BAMMD5\n","NA");
	return($BAMFile, $filteredDir);
}

#$snakemake .= "
#	input:
#		seqFile=$seqFile,
#		bismark_folder=$bismark_geneIndexDir
#	output:
#	
#";

