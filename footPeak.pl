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
use vars qw($opt_v $opt_g $opt_p $opt_d $opt_s $opt_k $opt_K $opt_n $opt_h $opt_t $opt_w $opt_L $opt_l $opt_o $opt_A $opt_G);
getopts("vg:p:d:s:k:K:n:ht:w:l:A:o:G:L:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
}

use myFootLib;
use FAlite;
use footPeakAddon;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
my $version = join("", @version);
if (defined $opt_v) {
	print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
	($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

##################
# 0. Check Sanity #
##################

my $date = getDate();
my $uuid = getuuid();
my $numThread = 1;
my ($footFolder) = $opt_n;
my ($usage) = check_sanity($footFolder);
my $logFile = "$footFolder/logFile.txt";

my ($opts, $outLog) = parse_footLoop_logFile($logFile, $date, $uuid, $footFolder, $version_small);
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
my $outDir 			= $opts->{n};
my $resDir        = $opts->{o};

#########################
# 1. Define Input/Names #
#########################

if (defined $opt_L) {
	print "label = $label\n";
}
elsif (-e "$footFolder/.LABEL") {
	($label) = `cat $footFolder/.LABEL`;
	chomp($label);
}
if (not defined $label or (defined $label and $label !~ /PCB/i)) {
	if ($resDir =~ /PCB\d+/i) {
		($label) = $resDir =~ /(PCB\d+(_BC\d+\_[A-Z]+)?)/i;
		$label = uc($label);
		$label =~ s/PCB0+(\d+)/PCB$1/ if $label =~ /PCB0+\d+/;
	}
	if (not defined $label and $resDir =~ /PCB[0\-_]*\d+/i) {
		($label) = $resDir =~ /(PCB[0\-_]*\d+(_BC\d+\_[A-Z]+)?)/i;
		$label = uc($label);
		$label =~ s/PCB[0\-_]+(\d+)/PCB$1/g if $label =~ /PCB[0\-_]+\d+/;
	}
	if (not defined $label and $label =~ /PCB/) {
		($label) = $resDir =~ /(PCB.{1,5}(_BC\d+\_[A-Z]+)?)\/?/;
	}
	if (not defined $label or (defined $label and $label !~ /PCB/i)) {
		die "Please make sure your output folder (-o) contain PCB(number) e.g. PCB12: 180202_PCB12_footpeak_output (no space/dash between PCB and number)\n\n";
	}
	else {
		LOG($outLog, "!label=$label\n");
	}
	if ($resDir =~ /_BC\d+_PFC\d+_\w+\./i) {
		my ($label2) = $resDir =~ /_(BC\w+)\./i;
		$label = "$label\_$label2";
		$label = uc($label);
		$label =~ s/PCB[0\-_]+(\d+)/PCB$1/g;
	}
}
system("echo $label > $resDir/.LABEL");
system("/bin/cp $seqFile $resDir");
LOG($outLog, "seqFile=$seqFile\n");
makedir("$resDir/.CALL") if not -d "$resDir/.CALL";

LOG($outLog, $usage . "faFile=$seqFile\nindexFile=$indexFile\norigFolder=$origFolder\noutdir=$outDir\n") and die unless defined $seqFile and defined $indexFile and defined $origFolder and -e $seqFile and -e $indexFile and -e $origFolder and defined $outDir;
#die $usage if mydefined($seqFile, $indexFile, $x, $y, $min, $groupsize, $minDis, $outDir) != 0;
if (not -d $outDir) {
	mkdir $outDir or LOG($outLog, "Cannot create directory (-n) $outDir: $!\n\n") and die;
}

# Get Real Coordinate of each Gene
my $SEQ = parse_indexFile_and_seqFile($indexFile, $seqFile, $outLog);

# get seq from faFile
# process peak

my @origFile = <$origFolder/*.orig>;
open (my $outLogAddon, ">", "$resDir/footLoop_addition_logFile.txt") or DIELOG($outLog, "Failed to write to $resDir/footLoop_addition_logFile.txt: $!\n");
close $outLogAddon;
my $out;
my $R;
my @total_Rscript = <$origFolder/*.R>;
my @total_png = <$origFolder/*.png>;
my %genes;
for (my $i = 0; $i < @origFile; $i++) {
	my ($peakFolder, $peakFilename) = getFilename($origFile[$i], "folderfull");
	$peakFilename =~ s/.orig$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	$genes{uc($gene)} ++;
}



LOG($outLog, "\n\n--------------$YW\n1. footPeak.pl processing $LGN" . scalar(keys %genes) . "$N genes (total file = $LGN" . scalar(@origFile) . "$N)\n\n");
for (my $i = 0; $i < @origFile; $i++) {
	my $Q = new Thread::Queue;
	my $peakFile = $origFile[$i];
	if (defined $opt_G and $peakFile !~ /$opt_G/i) {
		LOG($outLog, date() . "-> Skipped $LCY$peakFile$N as it doesn't contain $LGN-G $opt_G$N\n");
		next;
	}
#debug
#	next if $peakFile !~ /CALM3_Pos/;#CALM3_Pos/;#/(AIRN_PFC66_FORWARD_Neg|CALM3_Pos)/;
##
	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
	$peakFilename =~ s/.orig$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	DIELOG($outLog, "gene=$gene seq->gene{seq} isn't defined!\n") if not defined $SEQ->{$gene} or not defined $SEQ->{$gene}{seq};
	my ($totalPeak, $linecount, $peakcount, $total, %data, %end, %bad, %final, %group) = (0,0,0, scalar(@{$SEQ->{$gene}{seq}}));
	LOG($outLog, date() . "-> DEBUG FILENAME=$peakFilename, GENE=$gene\n","NA");
	DIELOG($outLog, "Died cannot get gene from peakfile $peakFile\n") unless defined $gene;
	DIELOG($outLog, "Cannot find sequence of gene=$LCY$gene$N in $seqFile!\n") if not defined $SEQ->{$gene};
	my %pk = ("CH" => 0, "CG" => 0, "GH" => 0, "GC" => 0);
	print "$LGN$i$N. peakFile=$LCY$peakFile$N, Gene = $LRD$gene$N total length = $LPR$total$N\n";
	foreach my $type (sort keys %pk) {
		my $outpeak = "PEAK$type";
		my $outnopk = "NOPK$type";
		open ($out->{$outpeak}, ">", "$resDir/.CALL/$peakFilename\_$window\_$threshold\_$type.PEAK") or DIELOG($outLog, "Failed to write to $resDir/.CALL/$peakFilename\_$window\_$threshold\_CG.PEAK: $!\n");
		open ($out->{$outnopk}, ">", "$resDir/.CALL/$peakFilename\_$window\_$threshold\_$type.NOPK") or DIELOG($outLog, "Failed to write to $resDir/.CALL/$peakFilename\_$window\_$threshold\_CG.NOPK: $!\n");
	}
	open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
	my ($l0, $t0) = (0,Benchmark->new());
	foreach my $type (sort keys %pk) {
		my $outpeak = "PEAK$type";
		my $outnopk = "NOPK$type";
		print {$out->{$outpeak}} "$peakFile\tPEAK\t$gene\t$type\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
		print {$out->{$outnopk}} "$peakFile\tNOPK\t$gene\t$type\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
	}

	LOG($outLog, date() . "$peakFilename\n");
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		next if $line =~ /^#/;
		my ($name, $type, $val) = split("\t", $line);
		$name = "footPeak.pl.name.UNDEF" if not defined $name or $name eq "" or $name =~ /^[ \t\n\s]*$/;
		$type = "footPeak.pl.type.UNDEF" if not defined $type or $type eq "" or $type =~ /^[ \t\n\s]*$/;
		$val = "footPeak.pl.val.UNDEF" if not defined $val or $val eq ""     or $val  =~ /^[ \t\n\s]*$/;
		if ($name eq "footPeak.pl.name.UNDEF" or $type eq "footPeak.pl.type.UNDEF" or $val eq "footPeak.pl.val.UNDEF") {
			LOG($outLog, "\n\n" . date() . "file=$YW$peakFilename$N line number $LGN$linecount$LRD ERROR$N: Skipped line error trying to parse$LCY tab delim$N :LGN\$name=$name $LPR\$type=$type $LCY\$val=$val$N, line:\n\n$line\n\n");
			next;
		}
		$val = [split("", $val)];
		$name = "$gene.$name";
		my ($count, $beg, $end) = (0,0,0);
		my $check = 0;
		if (defined $opt_G and $name !~ /$opt_G/i) {
			$check = 1;
			next;
		}
		my $comm = [$gene, $strand, $val, $name, $window, $threshold, $check];
		$Q->enqueue($comm);
		if ($Q->pending == 60) {
			LOG($outLog, date() . "-> Done $linecount\n") if $linecount % 300 == 0;
			$Q->end();
			$R = new Thread::Queue;
			my @threads; my $badthread = 0;
			for (my $thr = 0; $thr < $numThread; $thr++) {
				$threads[$thr] = threads->create(\&worker, $thr, $Q);
				if (defined $threads[$thr] and $threads[$thr] eq 1) {
					$badthread = 1;
					last;
				}
			}
			if ($badthread == 1) {
				for (my $thr = 0; $thr < $numThread; $thr++){
					$threads[$thr] = 1;
				}
			}
			for (my $thr = 0; $thr < $numThread; $thr++){
				$threads[$thr]->join();
			}
			$R->end();
			die if $badthread == 1;
			while ($R->pending) {
				my $result = $R->dequeue; next if not defined $result;
				my ($pkcurr, $outcurr) = @{$result};
				foreach my $type (sort keys %{$pkcurr}) {
					$pk{$type} += $pkcurr->{$type};
				}
				foreach my $type (sort keys %{$outcurr}) {
					print {$out->{$type}} $outcurr->{$type};
				}
				#last if $linecount > 100;
				#exit 0 if defined $check and $check == 1;
			}
			$Q = new Thread::Queue();
		}
		#last if $linecount > 100;
	}

	if (defined $Q->pending) {
		LOG($outLog, date() . "-> Done $linecount\n");
		$Q->end();
		$R = new Thread::Queue;
		my @threads; my $badthread = 0;
		for (my $thr = 0; $thr < $numThread; $thr++) {
			$threads[$thr] = threads->create(\&worker, $thr, $Q);
			if (defined $threads[$thr] and $threads[$thr] eq 1) {
				$badthread = 1;
				last;
			}
		}
		if ($badthread == 1) {
			for (my $thr = 0; $thr < $numThread; $thr++){
				$threads[$thr] = 1;
			}
		}
		for (my $thr = 0; $thr < $numThread; $thr++){
			$threads[$thr]->join();
		}
		$R->end();
		die if $badthread == 1;
		while ($R->pending) {
			my $result = $R->dequeue; next if not defined $result;
			my ($pkcurr, $outcurr) = @{$result};
			foreach my $type (sort keys %{$pkcurr}) {
				$pk{$type} += $pkcurr->{$type};
			}
			foreach my $type (sort keys %{$outcurr}) {
				print {$out->{$type}} $outcurr->{$type};
			}
			#last if $linecount > 100;
			#exit 0 if defined $check and $check == 1;
		}
	}
	my $to;
	foreach my $key (sort keys %pk) {
		next if $key =~ /NO$/;
		my $peak = defined $pk{$key} ? $pk{$key} : 0;
		$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
		LOG($outLog, "$key=$LGN$peak$N, ");
	}
	LOG($outLog, "\n\n" . date() . "\n\n");
	my $peakFilez = "$resDir/.CALL/$peakFilename\_$window\_$threshold\_CG.PEAK";
#   my ($input1, $faFile, $mygene, $minDis, $resDir, $minLen, $SEQ, $outLog) = @_;
	my $footPeakres = footPeakAddon::main(($peakFilez, $seqFile, $gene, $minDis, $resDir, $minLen, $SEQ, $version_small, $outLog));
	die "Failed footpeak\n" if defined $footPeakres and $footPeakres eq -1;
}
#system("/bin/cp $opts->{origDir}/footPeak_logFile.txt $resDir/");
#COMMENTCUT
#	my $t1 = Benchmark->new();
#	my ($td) = timestr(timediff($t1, $t0)) =~ /(\-?\d+\.?\d*) wallclock/;
#	my $persec = int(($linecount-$l0) / $td*10+0.5)/10;
#	LOG($outLog, "$peakFilename: Done $LGN$linecount$N in $LCY$persec$N lines per second; ");
#	my $to = 0;
#	foreach my $key (sort keys %pk) {
#		next if $key =~ /NO$/;
#		my $peak = defined $pk{$key} ? $pk{$key} : 0;
#		$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
#		LOG($outLog, "$key=$LGN$peak$N, ");
#	}
#	LOG($outLog, "total=$LCY$to$N\n");
##	print "PEAK: CH=$LGN$peakCH$N, CG=$LGN$peakCG$N, GH=$LCY$peakGH$N, GC=$LCY$peakGC$N\n";
##	print "NOPK: CH=$LGN$peakCH$N, CG=$LGN$peakCG$N, GH=$LCY$peakGH$N, GC=$LCY$peakGC$N\n";


###############
# SUBROUTINES #
###############

sub worker {
   my ($thread, $queue) = @_;
   my $tid = threads->tid;

   while ($queue->pending) {
		my $left = $queue->pending;
		return 0 if not defined $left;
#		print "Remaining: $left\n" if $left % 10 == 0;
      my $command = $queue->dequeue;
		my ($gene, $strand, $val, $name) = ($command->[0], $command->[1], $command->[2], $command->[3]);
		print "left; $left. tid=$tid doing $gene $name\n" if $left % 100 == 0;
		print "gene=$gene, strand=$strand, name=$name, Undefined seq!\n" and return 1 if not defined $SEQ;
		my ($peak, $peakCH, $peakCG, $peakGH, $peakGC) = getPeak(@{$command});
		#my ($gene, $strand, $nuc, $name, $out, $seq, $window, $threshold, $check) = @_;
		my ($CH, $CG, $GH, $GC) = ("","","","");
		for (my $i = 0; $i < @{$SEQ->{$gene}{seq}}; $i++) {
			$CH .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{CH}[$i] and $peak->{CH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CH}[$i] : "\t1";
			$CG .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{CG}[$i] and $peak->{CG}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CG}[$i] : "\t1";
			$GH .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{GH}[$i] and $peak->{GH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GH}[$i] : "\t1";
			$GC .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{GC}[$i] and $peak->{GC}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GC}[$i] : "\t1";
		}
		$CH = $peakCH != 0 ? "$name\tPEAK\t$gene\tCH\t$strand\t$CH\n" : "$name\tNOPK\t$gene\tCH\t$strand\t$CH\n";
		$CG = $peakCG != 0 ? "$name\tPEAK\t$gene\tCG\t$strand\t$CG\n" : "$name\tNOPK\t$gene\tCG\t$strand\t$CG\n";
		$GH = $peakGH != 0 ? "$name\tPEAK\t$gene\tGH\t$strand\t$GH\n" : "$name\tNOPK\t$gene\tGH\t$strand\t$GH\n";
		$GC = $peakGC != 0 ? "$name\tPEAK\t$gene\tGC\t$strand\t$GC\n" : "$name\tNOPK\t$gene\tGC\t$strand\t$GC\n";
		my ($pk, $out);
		if ($peakCH != 0) {$pk->{CH} ++; $out->{'PEAKCH'} = $CH;} else {$pk->{CHNO} ++; $out->{'NOPKCH'} = $CH;}
		if ($peakCG != 0) {$pk->{CG} ++; $out->{'PEAKCG'} = $CG;} else {$pk->{CGNO} ++; $out->{'NOPKCG'} = $CG;}
		if ($peakGH != 0) {$pk->{GH} ++; $out->{'PEAKGH'} = $GH;} else {$pk->{GHNO} ++; $out->{'NOPKGH'} = $GH;}
		if ($peakGC != 0) {$pk->{GC} ++; $out->{'PEAKGC'} = $GC;} else {$pk->{GCNO} ++; $out->{'NOPKGC'} = $GC;}
		$R->enqueue([$pk, $out]);
	}
	return 0;
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

sub getPeak {
	my ($gene, $strand, $nuc, $name, $window, $threshold, $check) = @_;
	#my $seq = $SEQ->{$gene};
   #CH=B(-)CDcdMN, CG=F(-)EeO
   #GH=J(-)GHghUV, GC=K(-)IiW
	#-_.TA
	#CDE123MNOPQR
	#GHI456UVWXYZ
	my $peak;
	my ($peakG, $peakC) = (0,0);
	my $con = makeCon();
	my ($peakCH, $peakCG, $peakGH, $peakGC) = (0,0,0,0);
	my ($chunkC, $chunkG); 
	@{$peak->{CH}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{CG}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{GH}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{GC}} = ("!") x (@{$SEQ->{$gene}{seq}});
	my $totalseqpos = (keys %{$SEQ->{$gene}{loc}{pos1}});
	my $totalseqneg = (keys %{$SEQ->{$gene}{loc}{neg1}});
	my $max = $totalseqpos > $totalseqneg ? $totalseqpos : $totalseqneg;
	# initialize
#	return($peak,0,0,0,0); #debug
	for (my $i = 0; $i < $max; $i++) {
		my ($jmin, $jmax) = $i == 0 ? (0, $window) : ($i, $i+1);
		for (my $j = $jmin; $j < $jmax; $j++) {
			if ($j < $totalseqpos - $window) {
				if ($i == 0) {
					my $nucpos1 = $nuc->[$SEQ->{$gene}{loc}{pos1}{$j}];
					$con = det_CG_con($con, $nucpos1, 1);
					LOG($outLog, "i=$i, j=$j: nuc.i=$LCY$nucpos1$N, nuc.j(nucpos1)=$nucpos1: file.orig value isn't one of the listed at $LGN" . __LINE__ . "$N\n") if $nucpos1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
					$peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CH");
#					$peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CG");
					$peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CG", $SEQ->{$gene}{loc}{pos1}{$j});
				}
				else {
					my ($nucpos0, $nucpos2) = ($nuc->[$SEQ->{$gene}{loc}{pos1}{$j-1}], $nuc->[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}]);
					$con = det_CG_con($con, $nucpos0, -1, $nucpos2, 1);
					$peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CH");
					$peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CG",$SEQ->{$gene}{loc}{pos1}{$j+$window-1});
				}
#				$chunkC .= $nucpos1;
			}
			if ($j < $totalseqneg - $window) {
				if ($i == 0) {
					my $nucneg1 = $nuc->[$SEQ->{$gene}{loc}{neg1}{$j}];
					$con = det_CG_con($con, $nucneg1, 1);
					LOG($outLog, "i=$i, j=$j: nuc.i=$LCY$nucneg1$N, nuc.j(nucneg1)=$nucneg1: file.orig value isn't one of the listed at $LGN" . __LINE__ . "$N\n") if $nucneg1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
					$peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j}] = det_peak($peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j}], $nucneg1, "GH");
					$peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j}] = det_peak($peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j}], $nucneg1, "GC");
				}
				else {
					my ($nucneg0, $nucneg2) = ($nuc->[$SEQ->{$gene}{loc}{neg1}{$j-1}], $nuc->[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}]);
					$con = det_CG_con($con, $nucneg0, -1, $nucneg2, 1);
					$peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}], $nucneg2, "GH");
					$peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}], $nucneg2, "GC");
				}
#				$chunkG .= $nucneg1;
			}
#=cut
		}
		if ($i == 0 and defined $check) {
#			print "i=$i, " . join("", @{$peak->{CH}}[0..$window-1]) . "\n";
		}
		my ($CG, $CH, $GC, $GH) = (0,0,0,0);
		($peak, $CG, $CH, $GC, $GH) = isPeak($i, $con, $peak, $window, $threshold, $SEQ->{$gene}, $totalseqpos, $totalseqneg, $check);
		$peakCH += $CH;
		$peakCG += $CG;
		$peakGC += $GC;
		$peakGH += $GH;
	}
	return ($peak, $peakCH, $peakCH, $peakGC, $peakGH);
}
sub det_peak { # time= 10 line/s per strand (49 l/s -> 39 l/s)
	my ($peak, $nuc, $type, $pos) = @_;
	$pos = "NA" if not defined $pos;
#	return 1; #debug
	return $peak if $peak ne '!';
	my @N =	$type eq 'CH' ? ('0', '0'   , 'cd', 'CDMNPQ', '12' , 'B\-' ) :
				$type eq 'CG' ? ('e', 'EQR' , 'cd', 'CDMNPQ', '123', 'BF\-') :
				$type eq 'GH' ? ('0', '0'   , 'gh', 'GHUVXY', '45' , 'J\-' ) :
				$type eq 'GC' ? ('i', 'IWZ' , 'gh', 'GHUVXY', '456', 'JK\-') : 
				die "det_peak TYPE IS NOT CH/CG/GH/GC? ($type) at line $LGN" . __LINE__ . "$N\n";
	return 7 if $nuc eq $N[0]; # conv CG
	return 5 if $nuc =~ /[$N[1]]/; # nonconverted CG
	return 6 if $nuc =~ /[$N[2]]/; #conv CH
	return 4 if $nuc =~ /[$N[3]]/; #nonconverted CH, or PQR = bad region plus MNOPQR = not C->C or C->T
	return 0 if $nuc =~ /[$N[4]]/; # 123 = C is located in bad region, and is C->C or C->T
	return 0 if $nuc =~ /[$N[5]]/; #C to deletion
	return ' ';
# 0 is bad or not data (6)
# 1 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak
}

sub det_CG_con { # time= 5 line/s (57 l/s -> 52 l/s)
	my ($con, @cons) = @_;
#	return $con; #debug
	for (my $i = 0; $i < @cons; $i+= 2) {
		my ($type, $val) = ($cons[$i], $cons[$i+1]);
		($con->{CG}{con} += $val and $con->{CG}{tot} += $val and next) if $type eq "e";
		($con->{GC}{con} += $val and $con->{GC}{tot} += $val and next) if $type eq "i";
		($con->{CH}{con} += $val and $con->{CH}{tot} += $val and next) if $type eq "c" or $type eq "d";
		($con->{GH}{con} += $val and $con->{GH}{tot} += $val and next) if $type eq "g" or $type eq "h";
#		($con->{CG}{tot} += $val and next) if $type eq "F" or $type eq "E" or $type eq "e" or $type eq "O" or $type eq "3" or $type eq "R";
		($con->{CG}{tot} += $val and next) if $type =~ /[FEeO3R]/;
		($con->{GC}{tot} += $val and next) if $type =~ /[KIiW6Z]/;
		($con->{CH}{tot} += $val and next) if $type =~ /[BCDcdMN12PQ]/;
		($con->{GH}{tot} += $val and next) if $type =~ /[JGHghUV45XY]/;
	}
	return $con;
}
sub makeCon {
	my $con;
	$con->{CH}{tot} = 0;
	$con->{CG}{tot} = 0;
	$con->{CH}{con} = 0;
	$con->{CG}{con} = 0;
	$con->{GH}{tot} = 0;
	$con->{GC}{tot} = 0;
	$con->{GH}{con} = 0;
	$con->{GC}{con} = 0;
	return $con;
}
sub isPeak {
	my ($i, $con, $peak, $window, $threshold, $seq, $totalseqpos, $totalseqneg, $check) = @_;
	my $CGcon = $con->{CG}{con} +$con->{CH}{con};
	my $CGtot = $con->{CG}{tot} +$con->{CH}{tot};
	my $CHcon = $con->{CH}{con};
	my $CHtot = $con->{CH}{tot};
	my $GCcon = $con->{GC}{con} +$con->{GH}{con};
	my $GCtot = $con->{GC}{tot} +$con->{GH}{tot};
	my $GHcon = $con->{GH}{con};
	my $GHtot = $con->{GH}{tot};
	my $peakCG = ($CGtot >= 5 and ($CGcon / $CGtot) > $threshold) ? 1 : 0;
	my $peakCH = ($CHtot >= 5 and ($CHcon / $CHtot) > $threshold) ? 1 : 0;
	my $peakGC = ($GCtot >= 5 and ($GCcon / $GCtot) > $threshold) ? 1 : 0;
	my $peakGH = ($GHtot >= 5 and ($GHcon / $GHtot) > $threshold) ? 1 : 0;
	if ($i < $totalseqpos - $window) {
		my $temp;
		my $min = $seq->{loc}{pos1}{$i};
		my $max = ($i < $totalseqpos - $window) ? $seq->{loc}{pos1}{$i + $window -1} : $seq->{loc}{pos1}{$totalseqpos-$window-1};
		if ($peakCH eq 1) {
			my $temp = join("", @{$peak->{CH}}[$min..$max]);
			$temp =~ tr/26/88/;
			if ($temp =~ /8$/) {
				$temp =~ s/8$/2/;
			}
			else {
				$temp =~ s/8([01345679A-Z\! ]+)$/2$1/i;# if $temp =~ /8([012345679A-Z\! ]+)$/i;
			}
			@{$peak->{CH}}[$min..$max] = split("", $temp);
			#die "$temp\n" if $peakCH eq 1;
	#		print "$i: $temp\n" if defined $check;
		}
		if ($peakCG eq 1) {
			my $temp = join("", @{$peak->{CG}}[$min..$max]);
			$temp =~ tr/2367/8989/;
			if ($temp =~ /8$/) {
				$temp =~ s/8$/2/;
			}
			elsif ($temp =~ /9$/) {
				$temp =~ s/9$/3/;
			}
			elsif ($temp =~ /8([014567A-Z\! \.]+)$/i) {
				$temp =~ s/8([014567A-Z\! ]+)$/2$1/i;
			}
			elsif ($temp =~ /9([014567A-Z\! \.]+)$/i) {
				$temp =~ s/9([01234567A-Z\! ]+)$/3$1/i;
			}
			@{$peak->{CG}}[$min..$max] = split("", $temp);
		}
#		print "i=$i, min=$min, max=$max, CHcon=$LGN$CHcon$N, CHtot=$LPR$CHtot$N, " . join("", @{$peak->{CH}}[$min..$max]) . "\n" if defined $check and $check eq 1;# and $peakCH eq 1;
	}
	if ($i < $totalseqneg - $window) {
		my $temp;
		my $min = $seq->{loc}{neg1}{$i};
		my $max = $i < $totalseqneg - $window ? $seq->{loc}{neg1}{$i + $window -1} : $seq->{loc}{neg1}{$totalseqneg-$window-1};
		if ($peakGH eq 1) {
			my $temp = join("", @{$peak->{GH}}[$min..$max]);
			$temp =~ tr/26/88/;
			if ($temp =~ /8$/) {
				$temp =~ s/8$/2/;
			}
			else {
				$temp =~ s/8([012345679A-Z\! \.]+)$/2$1/i;# if $temp =~ /8([012345679A-Z\! ]+)$/i;
			}
			@{$peak->{GH}}[$min..$max] = split("", $temp);
		}
		if ($peakGC eq 1) {
			my $temp = join("", @{$peak->{GC}}[$min..$max]);
			$temp =~ tr/2367/8989/;
			if ($temp =~ /8$/) {
				$temp =~ s/8$/2/;
			}
			elsif ($temp =~ /9$/) {
				$temp =~ s/9$/3/;
			}
			elsif ($temp =~ /8([01234567A-Z\! \.]+)$/i) {
				$temp =~ s/8([01234567A-Z\! ]+)$/2$1/i;
			}
			elsif ($temp =~ /9([01234567A-Z\! \.]+)$/i) {
				$temp =~ s/9([01234567A-Z\! ]+)$/3$1/i;
			}
			@{$peak->{GC}}[$min..$max] = split("", $temp);
		}
	}
#	return ($peak, 0,0,0,0);
	return ($peak, $peakCG, $peakCH, $peakGC, $peakGH);
}
sub mydefined {
	my @arr = @_;
	for (my $i = 0; $i < @arr; $i++) {
		(print "$i: mydefined failed\n" and return 1) if not defined $arr[$i];
	}
	return 0;
}
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
footLoop samFile     : $defOpts->{samFile}
footLoop samFixed    : $defOpts->{samFixed}
footLoop samFixedMD5 : $defOpts->{samFixedMD5}
>Options from footLoop.pl logfile=$logFile\n";
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
	my ($logFile, $date, $uuid, $footFolder, $version_small) = @_;
	my @line = `cat $logFile`;
	my $paramsFile = "$footFolder/.PARAMS";
	my ($defOpts, $usrOpts, $usrOpts2) = set_default_opts();
	my $inputFolder = $defOpts->{n}; 

	my @parline = `cat $paramsFile`;
	foreach my $parline (@parline) {
		if ($parline =~ /footLoop.pl,geneIndexFile,/) {
			($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,geneIndexFile,/) {
			($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,seqFile,/) {
			($defOpts->{seqFile}) = $parline =~ /seqFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFile,/) {
			($defOpts->{samFile}) = $parline =~ /samFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFixedFile,/) {
			($defOpts->{samFixed}) = $parline =~ /samFixedFile,(.+)$/;
		}
		if ($parline =~ /footLoop.pl,samFixedFileMD5,/) {
			($defOpts->{samFixedMD5}) = $parline =~ /samFixedFileMD5,(.+)$/;
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
			$value = $footFolder . "/$value";
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
#			$defOpts->{geneIndexFile} = $footFolder . "/" .  getFilename($defOpts->{geneIndexFile});
#		}
		if ($line =~ /^!\w+=/) {
			my ($param, $value) = $line =~ /^!(\w+)=(.+)$/;
			my $param2 = defined $param ? $param : "__UNDEF__";
			my $value2 = defined $value ? $value : "__UNDEF__";
			if ($value =~ /\//) {
				if ($value =~ /\/?\.0_orig\w{32}/) {
					($value) = $value =~ /^.+(\/?\.0_orig_\w{32})/; 
					$value = $footFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
				if ($value =~ /\/\.geneIndex/) {
					($value) = $value =~ /^.+(\/\.geneIndex.+)$/; 
					$value = $footFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
				else {
					($value) = getFilename($value, 'full');
					$value = $footFolder . "/$value";
					die "Died at line=line, param=$param, value=?\n" if not defined $value;
				}
			}
			print "$param = $value\n";
			print "Cannot parse param=$param2 and value=$value2 from line=$line\n" and die if not defined $param or not defined $value;
			print "$param file $value does not exist!\n" and die if $value =~ /\/+/ and not -e $value;
			($defOpts->{$param}) = $value if $param ne "n";
		}
	}
	$defOpts->{o} = $defOpts->{n} if not defined $opt_o;
	makedir($defOpts->{o}) if not -d $defOpts->{o};
	#die "opt = $opt_o = $defOpts->{o}\n";
	open ($outLog, ">", "$defOpts->{o}/footPeak_logFile.txt") or print "Failed to write to $defOpts->{o}/footPeak_logFile.txt: $!\n" and die;
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
			print "parse_runscript: opt $opt in footLoop logFile.txt Run Script does not look like a getopt options (not a '-a')!\n" and next;
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
		'y' => '0'       ,'K' => '2'       ,'L' => ''       ,'A' => ''        ,
		'o' => 'RESULTS' ,'G' => ''
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
	LOG($outLog, "Index File = $indexFile\n");
	foreach my $line (@indexLine) {
		chomp($line);
		my ($chr, $beg, $end, $def, $zero, $strand) = split("\t", $line);
		$def = uc($def);
		$SEQ->{$def}{coor} = "$chr\t$beg\t$end\t$def\t$zero\t$strand";
#		LOG($outLog, "indexFile:\tSEQ -> {gene=$def}{coor} is defined!\n");
		LOG($outLog, "def=$def, coor=$chr, $beg, $end, $def, $zero, $strand\n");
	}

	open (my $in, "<", $seqFile) or DIELOG($outLog, "Failed to read from $LCY$seqFile$N: $!\n");
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def; $def =~ s/^>//; $def = uc($def);
		my @seq = split("", $entry->seq);
		$SEQ->{$def}{seq} = \@seq;
		$SEQ->{$def}{loc} = findCGPos(\@seq);
		LOG($outLog, "seqFile:\tSEQ -> {gene=$def}{seq} and {loc} is defined!\n");
	}
	close $in;
	return $SEQ;
}

sub check_sanity {
	my ($footFolder) = @_;
	my ($usage, $usage_long) = usage();
	if (defined $opt_h or not defined $footFolder or not -d $footFolder) {
		print $usage;
		print $usage_long if defined $opt_h;
		print "${LRD}ERROR:$N Please define -n (output dir from footLoop.pl)\n" if not defined $footFolder or not -d $footFolder;
		print "${LRD}ERROR:$N Please define -o (output dir of footPeak.pl)\n" if not defined $opt_o;
		die "${YW}----------------------$N\n\n";
	}
	LOG($outLog, "Please run footloop.pl first! ($footFolder/logFile.txt does not exists)\n") and die if not -e "$footFolder/logFile.txt";

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
