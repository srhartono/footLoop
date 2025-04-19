#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";

}


use myFootLib;
use FAlite;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);

my ($thisfileName) = getFilename($0, "fullname");
my ($thismd5) = getMD5_simple($0);

my ($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $version_small, $debug) = @ARGV;
die "\nUsage: $YW$0$N ${LCY}indexFile$N ${LGN}seqFile$N ${LCY}i$N ${LGN}origFile$N ${LCY}outDir$N ${LGN}window$N ${LCY}threshold$N ${LGN}totalorigFile$N ${LCY}label$N ${LGN}genewant$N ${LCY}minDis$N ${LGN}minLen$N ${LCY}version_small$N\n\n" unless @ARGV >= 13;
undef $genewant if defined($genewant) and $genewant eq -1;

my ($origFolder, $origFilename) = getFilename($origFile,"folderfull");

my $outLogFile = "$outDir/.footPeak_sbatch/$origFilename.footPeak_sbatch.log.txt";
open (my $outLog, ">", $outLogFile) or die "Failed to write to $LCY$outLogFile$N: $!\n";

my $SEQ = parse_indexFile_and_seqFile($indexFile, $seqFile, $outLog);

if (defined $debug) {
	LOG($outLog, "\n" . date() . "${LGN}SUCCESS!!$N Dry run $0 ran successfully\n\n");
	exit 0;
}
run_footPeak($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog);

LOG($outLog, date() . "Done!\n\n$outLogFile\n\n");
close $outLog;

###############
# SUBROUTINES #
###############

sub run_footPeak {
	my ($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog) = @_;
	my $resDir = $outDir;
	my $out;
	my $peakFile = $origFile;
	LOG($outLog, date() . "$LGN$origFileInd$N/$LGN$totalorigFile$N:$LPR $0 processing$N $YW$peakFile$N\n");

	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
	$peakFilename =~ s/.filtered.gz$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);

	
	DIELOG($outLog, "gene=$gene seq->gene{seq} isn't defined!\n") if not defined $SEQ->{$gene} or (defined $SEQ->{$gene} and not defined $SEQ->{$gene}{seq});

	my ($totalPeak, $linecount, $peakcount, $total, %data, %end, %bad, %final, %group) = (0,0,0, scalar(@{$SEQ->{$gene}{seq}}));
	
	DIELOG($outLog, "Died cannot get gene from peakfile $peakFile\n") unless defined $gene;
	DIELOG($outLog, "Cannot find sequence of gene=$LCY$gene$N in $seqFile!\n") if not defined $SEQ->{$gene};
	my %pk = ("CH" => 0, "CG" => 0, "GH" => 0, "GC" => 0);
	my @seq = @{$SEQ->{$gene}{seq}};

	foreach my $dinuc (sort keys %pk) {
		open (my $outpeak, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK: $!\n");
		open (my $outpeakout, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.out") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK: $!\n");
		open (my $outnopk, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK: $!\n");
		open (my $outnopkout, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.out") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK: $!\n");
		print $outpeak "$peakFile\tPEAK\t$gene\t$dinuc\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
		print $outnopk "$peakFile\tNOPK\t$gene\t$dinuc\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
		close $outpeak;
		close $outpeakout;
		close $outnopk;
		close $outnopkout;

		open (my $outpeakconvperc, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.convperc") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.convperc: $!\n");
		close $outpeakconvperc;
		open (my $outnopkconvperc, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.convperc") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.convperc: $!\n");
		close $outnopkconvperc;


		my $outPEAK_GENOME = "$resDir/PEAKS_GENOME/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.genome.bed";
		my $outPEAK_LOCAL = "$resDir/PEAKS_LOCAL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.local.bed";
		open (my $outPEAK_GENOMEout, ">", $outPEAK_GENOME) or LOG($outLog, "\tFailed to write into $outPEAK_GENOME: $!\n")  and exit 1;
		open (my $outPEAK_LOCALout, ">", $outPEAK_LOCAL) or LOG($outLog, "\tFailed to write into $outPEAK_LOCAL: $!\n")  and exit 1;
		close $outPEAK_GENOMEout;
		close $outPEAK_LOCALout;
		print "PEAKFILE: $LPR$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK$N\n"
	}
	my ($peakFilesize) = -s $peakFile; 
	LOG($outLog, "$LGN$origFileInd$N. size=$LGN$peakFilesize$N Gene = $LRD$gene$N total length = $LPR$total$N\n","NA") if $peakFilesize <= 20;

	my ($peakFilelinecount) = `zcat $peakFile|wc -l` =~ /^(\d+)$/;
	LOG($outLog, "$LGN$origFileInd$N/$LGN" . $totalorigFile .  "$N $LGN$peakFilename$N $LGN$peakFilelinecount$N size=$LGN$peakFilesize$N Gene = $LRD$gene$N total length = $LPR$total$N totalread=$LCY$peakFilelinecount$N\n") if $peakFilelinecount >= 5;
	my $in1;
	if ($peakFile =~ /.gz$/) {
		open ($in1, "zcat $peakFile|") or die "Cannot read from $peakFile: $!\n";
	}
	else {
		open ($in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
	}

	while (my $line = <$in1>) {
		chomp($line);
		LOG($outLog, date() . "${LPR}$0$N: Done $LGN$linecount$N\n") if $linecount % 100 == 0;
		$linecount ++;
		next if $line =~ /^#/;
		my ($name, $type, $val) = split("\t", $line);
		$name = "$0.name.UNDEF" if not defined $name or $name eq "" or $name =~ /^[ \t\n\s]*$/;
		$type = "$0.type.UNDEF" if not defined $type or $type eq "" or $type =~ /^[ \t\n\s]*$/;
		$val = "$0.val.UNDEF" if not defined $val or $val eq ""     or $val  =~ /^[ \t\n\s]*$/;
		if ($name eq "$0.name.UNDEF" or $type eq "$0.type.UNDEF" or $val eq "$0.val.UNDEF") {
			LOG($outLog, "\n\n" . date() . "file=$YW$peakFilename$N line number $LGN$linecount$LRD ERROR$N: Skipped line error trying to parse$LCY tab delim$N :LGN\$name=$name $LPR\$type=$type $LCY\$val=$val$N, line:\n\n$line\n\n","NA");
			next;
		}
		$name = "$gene.$name";
		my ($count, $beg, $end) = (0,0,0);
		my $check = 0;

		if (defined $genewant and $name !~ /$genewant/i) {$check = 1; next;} 
		my @coor = split("\t", $SEQ->{$gene}{coor});
		my ($chr0, $beg0, $end0, $name0, $val0, $strand0) = @coor;

		my @refs = @{$SEQ->{$gene}{seq}};
		my $refs = join("", @refs);
		my $Ccon = $val =~ tr/cde/cde/;
		my $Ctot = $refs =~ tr/C/C/;
		my $Gcon = $val =~ tr/ghi/ghi/;
		my $Gtot = $refs =~ tr/G/G/;
		my $CTperc = $Ctot == 0 ? 0 : int(1000*$Ccon / $Ctot+0.5)/10;

		my $GAperc = $Gtot == 0 ? 0 : int(1000*$Gcon / $Gtot+0.5)/10;
		$val = [split("", $val)];

		my ($pkcurr) = dothis($gene, $strand, $val, $name, $window, $threshold, $minLen, $check);

		foreach my $dinuc (sort keys %{$pkcurr}) {
			$pk{$dinuc} += $pkcurr->{$dinuc}{PEAK}{TOTAL};
			if ($pkcurr->{$dinuc}{PEAK}{TOTAL} > 0) {
				open (my $outpeak, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK: $!\n");
				print $outpeak "$name\t$pkcurr->{$dinuc}{PEAK}{CONV}\n";
				close $outpeak;

				open (my $outpeakconvperc, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.convperc") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.convperc: $!\n");
				print $outpeakconvperc "$name\t$pkcurr->{$dinuc}{PEAK}{CONVPERC}\n";
				close $outpeakconvperc;

				open (my $outpeakout, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.out") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK: $!\n");
				print $outpeakout "$name\t$pkcurr->{$dinuc}{PEAK}{CONV}\n";
				close $outpeakout;
			}
			if ($pkcurr->{$dinuc}{NOPK}{TOTAL} > 0) {
				open (my $outnopk, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK: $!\n");
				print $outnopk "$name\t$pkcurr->{$dinuc}{NOPK}{CONV}\n";
				close $outnopk;

				open (my $outnopkconvperc, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.convperc") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.convperc: $!\n");
				print $outnopkconvperc "$name\t$pkcurr->{$dinuc}{NOPK}{CONVPERC}\n";
				close $outnopkconvperc;

				open (my $outnopkout, ">>", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.out") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK: $!\n");
				print $outnopkout "$name\t$pkcurr->{$dinuc}{NOPK}{CONV}\n";
				close $outnopkout;
			}
			my ($junk, $readname) = $name =~ /^(\w+)\.(.+)$/;
			if (defined $pkcurr->{$dinuc}{PEAK}{COOR}) {
				my $outPEAK_GENOME = "$resDir/PEAKS_GENOME/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.genome.bed";
				my $outPEAK_LOCAL = "$resDir/PEAKS_LOCAL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.local.bed";
				open (my $outPEAK_GENOMEout, ">>", $outPEAK_GENOME) or LOG($outLog, "\tFailed to write into $outPEAK_GENOME: $!\n")  and exit 1;
				open (my $outPEAK_LOCALout, ">>", $outPEAK_LOCAL) or LOG($outLog, "\tFailed to write into $outPEAK_LOCAL: $!\n")  and exit 1;
				foreach my $beg (sort{$a <=> $b} keys %{$pkcurr->{$dinuc}{PEAK}{COOR}}) {
					my $end = $pkcurr->{$dinuc}{PEAK}{COOR}{$beg};
					my $conv = $pkcurr->{$dinuc}{PEAK}{CONVPERCBED}{$beg};
					my $end1 = $end + $beg0;
					my $beg1 = $beg + $beg0;
					my $file = "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK";
					print $outPEAK_GENOMEout "$chr0\t$beg1\t$end1\t$name\t$conv\t$strand0\t$file\n";
					print $outPEAK_LOCALout "$name\t$beg\t$end\t$conv\n";
				}
				close $outPEAK_GENOMEout;
				close $outPEAK_LOCALout;
			}
		}
		if (defined $genewant) {
			last if $name =~ /$genewant/i;
		}
	}
	LOG($outLog, date() . "${LPR}$0$N: Parsed $LGN$linecount$N lines\n\n");
	my $to;
	foreach my $key (sort keys %pk) {
		next if $key =~ /NO$/;
		my $peak = defined $pk{$key} ? $pk{$key} : 0;
		$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
		LOG($outLog, "$key=$LGN$peak$N, ","NA");
	}
	LOG($outLog, "\n\n" . date() . "\n\n","NA");
	foreach my $dinuc (sort keys %pk) {
		my $peakoutFile = "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.out";
		my $nopkoutFile = "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.out";
		my ($totalpeak) = (-e $peakoutFile) ? `wc -l $peakoutFile` : " 0 $peakoutFile"; chomp($totalpeak);
		my ($totalnopk) = (-e $nopkoutFile) ? `wc -l $nopkoutFile` : " 0 $nopkoutFile"; chomp($totalnopk);
		open (my $outpeakoutwcl, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.wcl") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.PEAK.wcl: $!\n");
		print $outpeakoutwcl "$totalpeak\n";
		close $outpeakoutwcl;
		open (my $outnopkoutwcl, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.wcl") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_$dinuc.NOPK.wcl: $!\n");
		print $outnopkoutwcl "$totalnopk\n";
		close $outnopkoutwcl;
	}
}

sub make_total_hash {
	my $total;
	my @types = qw(CH CG GH GC);
	foreach my $rconvType (@types) {
		$total->{$rconvType}{peak} = 0; 
		$total->{$rconvType}{nopk} = 0;
		$total->{$rconvType}{total} = 0;
	}
	return($total);
}

sub parse_peak {
	my ($ARG, $bad, $minDis, $minLen, $outLog) = @_;
	my ($name, $isPeak, $mygene, $rconvType, $readStrand, @val) = split("\t", $ARG);
	my %bad = %{$bad} if defined $bad;
	my $name_want = "DEBUG";
	shift(@val) if $val[0] eq "";
	my $peaks;
	my %peak; $peak{curr} = 0; 
	my $Length = @val; 
	my $print = "name=$name, isPeak = $isPeak, Total length = $Length\n";
	my ($edge1) = join("", @val) =~ /^(0+)[\.1-9A-Za-z]/;
	$edge1 = defined $edge1 ? length($edge1) : 0;
	my ($edge2) = join("", @val) =~ /[\.1-9A-Za-z](0+)$/;
	$edge2 = defined $edge2 ? @val-length($edge2) : @val;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		if ($i % 100 == 0) {$print .= "\n$YW" . $i . "$N:\t";}
		if ($val[$i] =~ /[89]/) {
			$peak{beg} = $i if $peak{curr} == 0;
			$print .= "${LPR}$val[$i]$N" if $peak{curr} == 0;
			$print .= "${LRD}$val[$i]$N" if $peak{curr} == 1;
			$peak{curr} = 1;
		}
		elsif ($val[$i] =~ /[23]/) {
			$peak{end} = $i+1;
			push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
			undef $peak{beg}; undef $peak{end};
			$peak{curr} = 0;
			$val[$i] =~ tr/23/89/;
			$print .= "${LPR}$val$N";
		}
		else {
			$print .= "EDGE1" if $i == $edge1;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[57]$/;
			$print .= "." if $val[$i] eq 1;
			$print .= "x" if $val[$i] eq 0;
			$print .= "EDGE2" if $i == $edge2 - 1;
		}
	}
	my (%nopk, @peak);
	$readStrand = $rconvType =~ /^C/ ? 0 : $rconvType =~ /^G/ ? 16 : 255;
	my %peak2;
	$print .= "\n";
	print "$print" if $name_want eq $name;
	if (defined $peak{peak}) {
		foreach my $peak (sort @{$peak{peak}}) {
			my ($beg, $end) = split("-", $peak);
			my $checkBad = 0;
			foreach my $begBad (sort keys %{$bad->{$readStrand}}) {
				my $endBad = $bad->{$readStrand}{$begBad};
				if ($beg >= $begBad and $beg <= $endBad and $end >= $begBad and $end <= $endBad) {
					LOG($outLog, date() . "\t\t$LGN YES$N peak=$beg-$end, bad=$begBad-$endBad\n","NA") if $isPeak eq "PEAK";
					$checkBad = 1; last;
				}
			}
			if ($checkBad != 1) {
				next if not defined $bad->{$readStrand};
				foreach my $begBad (sort keys %{$bad->{$readStrand}}) {
					my $endBad = $bad->{$readStrand}{$begBad};
						my @valz = @val;
						my ($goodC, $badC) = (0,0);
						for (my $m = $beg; $m < $end; $m++) {
							if ($m >= $begBad and $m <= $endBad) {
								$badC ++ if $valz[$m] =~ /[2389]/;
							}
							else {
								$goodC ++ if $valz[$m] =~ /[2389]/;
							}
						}
						if ($goodC < 5 and $badC >= 9) {
							$checkBad = 1; last;
						}
				}
			}
			
			if ($checkBad == 1) {
				$print .= "\tCheckBad; Peak Not: $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
			elsif ($end - $beg < $minLen) {
				$print .= "\tend-$end < $minLen; Peak Not: $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
			elsif ($beg < $edge2 - 100 and $end > 100 + $edge1) {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Used: $LGN$peak$N\n";
				push(@peak, "$beg-$end");
				push(@{$peak2{peak}}, $peak);
			}
			else {
				$print .= "\tbeg=$beg, end=4end, peak Not: $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
		}
	}
	my $totalpeak = scalar(@peak);
	my @val2 = @val;
	if ($totalpeak > 0) {
		for (my $i = 0; $i < @val; $i++) {
			my $val = $val[$i];
			$val2[$i] = $val;
			if ($val =~ /^(8|9)$/ and defined $nopk{$i}) { 
				$val2[$i] = 7 if $val eq 9;
				$val2[$i] = 6 if $val eq 8;
			}
		}
	}
	$print .= "$name\t$totalpeak\n" if $isPeak eq "PEAK";
	@val = @val2;
	$print .= "\n\nVAL2: Total length = $Length\n";
	($edge1) = join("", @val) =~ /^(0+)[\.1-9A-Za-z]/;
	$edge1 = defined $edge1 ? length($edge1) : 0;
	($edge2) = join("", @val) =~ /[\.1-9A-Za-z](0+)$/;
	$edge2 = defined $edge2 ? @val-length($edge2) : @val;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		if ($i % 100 == 0) {$print .= "\n$YW" . $i . "$N:\t";}
		if ($val[$i] =~ /[89]/) {
			$peak{beg} = $i if $peak{curr} == 0;
			$print .= "${LPR}$val[$i]$N" if $peak{curr} == 0;
			$print .= "${LRD}$val[$i]$N" if $peak{curr} == 1;
			$peak{curr} = 1;
		}
		elsif ($val[$i] =~ /[23]/) {
			$peak{end} = $i+1;
			push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
			undef $peak{beg}; undef $peak{end};
			$peak{curr} = 0;
			$val[$i] =~ tr/23/89/;
			$print .= "${LPR}$val$N";
		}
		else {
			$print .= "EDGE1" if $i == $edge1;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[57]$/;
			$print .= "." if $val[$i] eq 1;
			$print .= "x" if $val[$i] eq 0;
			$print .= "EDGE2" if $i == $edge2 - 1;
		}
	}
	$print .= "\n";
	print "$print" if $name_want eq $name;

	return ($name, \@val2, $totalpeak, $peak2{peak});
}

sub find_lots_of_C {
	my ($seqFile, $mygene, $outLog) = @_;
	my %seq;
	LOG($outLog, "\n------------\n${YW}2. Getting polyC (or G) regions from sequence file $LCY$seqFile$N\n\n\n","NA");
	open(my $SEQIN, "<", $seqFile) or LOG($outLog, date() . "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!") and exit 1;
	my $fasta = new FAlite($SEQIN);
	my %lotsOfC;

	
	while (my $entry = $fasta->nextEntry()) {
	   my $gene = uc($entry->def); $gene =~ s/^>//;
	   my $seqz = uc($entry->seq);
		next if $gene ne $mygene;
	   LOG($outLog, date . "gene=$LGN$gene$N\n","NA");
	   
	
		my $minlen = 6;
	   my $seqz2 = $seqz;
	   while ($seqz2 =~ /(C){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_C) = length($curr);
	      my ($next_C) = $next =~ /^(C+)[AGTN]*$/;
	      $next_C = defined $next_C ? length($next_C) : 0;
	      my ($beg_C) = defined $prev ? length($prev) : 0;
	      my ($end_C) = $curr_C + $next_C + $beg_C;
	      my $length = $curr_C + $next_C;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_C to $end_C ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n","NA","NA");
	      $lotsOfC{$gene} .= "C;$beg_C;$end_C,";
	   }
		$seqz2 = "";
	   $seqz2 =$seqz;
	   while ($seqz2 =~ /(G){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_G) = length($curr);
	      my ($next_G) = $next =~ /^(G+)[ACTN]*$/;
	      $next_G = defined $next_G ? length($next_G) : 0;
	      my ($beg_G) = defined $prev ? length($prev) : 0;
	      my ($end_G) = $curr_G + $next_G + $beg_G;
	      my $length = $curr_G + $next_G;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_G to $end_G ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n","NA","NA");
	      $lotsOfC{$gene} .= "G;$beg_G;$end_G,";
	   }
		my $bad;
		if (defined $lotsOfC{$gene}) {
			LOG($outLog, date() . "Region of Bad C's:\n","NA");
			open (my $outBadC, ">", "$seqFile.$mygene.badc.bed") or DIELOG($outLog, "Failed to write to $seqFile.$mygene.badc.bed: $!\n");
			$lotsOfC{$gene} =~ s/,$//;
			my @lotsOfC = split(",", $lotsOfC{$gene});
			foreach my $coor (@lotsOfC) {
				my ($nuc, $beg, $end) = split(";", $coor);
				my $strands = $nuc eq "C" ? 0 : $nuc eq "G" ? 16 : 255;
				my $strandsymbol = $strands eq 0 ? "+" : "-";
				$bad->{$strands}{$beg} = $end-1;
				LOG($outLog, date() . "$strands: coor=$coor, nuc=$nuc, beg=$beg, end=$end, beg-end=$beg-$bad->{$strands}{$beg}\n","NA");
				print $outBadC "$gene\t$beg\t$end\t$coor\t" . ($end-$beg) . "\t$strandsymbol\n";
			}
			LOG($outLog, "\n${LGN}SUCCESS!$N PolyC (or G) bed file:\n$LCY$seqFile.$mygene.badc.bed$N\n\n","NA");
			return $bad;
		}
	}
	return;
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

sub dothis {
   my @command = @_;
   my ($gene, $strand, $val, $name, $window, $threshold, $minLen, $check) = @command;
   print "gene=$gene, strand=$strand, name=$name, Undefined seq!\n" and return 1 if not defined $SEQ;
   my ($peak) = getPeak2(@command);
   my ($CH, $CG, $GH, $GC) = ("","","","");
	my $peakGCz;
	my @dinuc = qw(CH CG GH GC);
  	for (my $i = 0; $i < @{$SEQ->{$gene}{seq}}; $i++) {
		$peakGCz .= ($peak->{GC}{conv}[$i] !~ /^[\! ]$/) ? $peak->{GC}{conv}[$i] : 1;
	}
	foreach my $beg (sort keys %{$peak->{GC}{peakcoor}}) {
		my $end = $peak->{GC}{peakcoor}{$beg};
	}
	foreach my $beg (sort keys %{$peak->{GC}{nopkcoor}}) {
		my $end = $peak->{GC}{nopkcoor}{$beg};
	}
   my ($pk, $out);
	foreach my $dinuc (@dinuc[0..@dinuc-1]) {
		$pk->{$dinuc}{PEAK}{TOTAL} = 0;
		$pk->{$dinuc}{NOPK}{TOTAL} = 0;
   	for (my $i = 0; $i < @{$SEQ->{$gene}{seq}}; $i++) {
	  	   die "Undefined i=$LGN$i$N \$val->[\$i]\n" if not defined $val->[$i];
			$peak->{$dinuc}{convfinal} .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{$dinuc}{conv}[$i] and $peak->{$dinuc}{conv}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{$dinuc}{conv}[$i] : "\t1";
		}
		$peak->{$dinuc}{convfinal} =~ s/^\t//;
	   $pk->{$dinuc}{convfinal} = $peak->{$dinuc}{tot} != 0 ? "$name\tPEAK\t$gene\t$dinuc\t$strand\t$peak->{$dinuc}{convfinal}\n" : "$name\tNOPK\t$gene\t$dinuc\t$strand\t$peak->{$dinuc}{convfinal}\n";
	   if ($peak->{$dinuc}{tot} != 0) {
			$pk->{$dinuc}{PEAK}{TOTAL} ++;
			$pk->{$dinuc}{PEAK}{CONVPERCBED} = $peak->{$dinuc}{convpercbed} if defined $peak->{$dinuc}{convpercbed};
			$pk->{$dinuc}{PEAK}{CONV} = $peak->{$dinuc}{convfinal};
			$pk->{$dinuc}{PEAK}{COOR} = $peak->{$dinuc}{peakcoor} if defined $peak->{$dinuc}{peakcoor};
			$pk->{$dinuc}{PEAK}{CONVPERC} = $peak->{$dinuc}{convperc};
		}
		else {	
			$pk->{$dinuc}{NOPK}{TOTAL} ++;
			$pk->{$dinuc}{NOPK}{CONV} = $peak->{$dinuc}{convfinal};
			$pk->{$dinuc}{NOPK}{COOR} = $peak->{$dinuc}{nopkcoor} if defined $peak->{$dinuc}{nopkcoor};
			$pk->{$dinuc}{NOPK}{CONVPERC} = $peak->{$dinuc}{convperc};
		}
	}
   return($pk);
}

sub getPeak2 {
   my ($gene, $strand, $nuc, $name, $window, $threshold, $minLen, $check) = @_;
	my $con = makeCon();
	my ($peakG, $peakC, $peakCH, $peakCG, $peakGH, $peakGC) = (0,0,0,0,0,0);
	my ($peak, $chunkC, $chunkG);
	@{$peak->{CH}{conv}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{CG}{conv}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{GH}{conv}} = ("!") x (@{$SEQ->{$gene}{seq}});
	@{$peak->{GC}{conv}} = ("!") x (@{$SEQ->{$gene}{seq}});
	my $totalseqpos = (keys %{$SEQ->{$gene}{loc}{pos1}});
	my $totalseqneg = (keys %{$SEQ->{$gene}{loc}{neg1}});
	my $max = $totalseqpos > $totalseqneg ? $totalseqpos : $totalseqneg;
	my @vals = @{$nuc};
	my @refs = @{$SEQ->{$gene}{seq}};
	my %CG;
	my %pos;
	my $CGind = 0;
	my $CHind = 0;
	my $GCind = 0;
	my $GHind = 0;
	$pos{C}{CGN}[0] = 0;
	$pos{C}{NCG}[0] = 0;
	$pos{G}{GCN}[0] = 0;
	$pos{G}{NGC}[0] = 0;
	$pos{C}{CHN}[0] = 0;
	$pos{C}{NCH}[0] = 0;
	$pos{G}{GHN}[0] = 0;
	$pos{G}{NGH}[0] = 0;
	@{$pos{C}{CGseq}} = ();
	@{$pos{G}{GCseq}} = ();
	@{$pos{C}{CHseq}} = ();
	@{$pos{G}{GHseq}} = ();
	$pos{C}{CGNseq}[0] = "";
	$pos{C}{NCGseq}[0] = "";
	$pos{G}{GCNseq}[0] = "";
	$pos{G}{NGCseq}[0] = "";
	$pos{C}{CHNseq}[0] = "";
	$pos{C}{NCHseq}[0] = "";
	$pos{G}{GHNseq}[0] = "";
	$pos{G}{NGHseq}[0] = "";
	for (my $i = 0; $i < @refs; $i++) {
		my $ref0 = $i == 0 ? "N" : $refs[$i-1];
		my $ref1 = $refs[$i+0];
		my $ref2 = $i == @refs - 1 ? "N" : $refs[$i+1];
		my $val1 = $vals[$i+0];
	 	($peak->{CH}{conv}[$i]) = det_peak2($peak->{CH}{conv}[$i], $val1, "CH");
	 	($peak->{CG}{conv}[$i]) = det_peak2($peak->{CG}{conv}[$i], $val1, "CG");
	 	($peak->{GH}{conv}[$i]) = det_peak2($peak->{GH}{conv}[$i], $val1, "GH");
	 	($peak->{GC}{conv}[$i]) = det_peak2($peak->{GC}{conv}[$i], $val1, "GC");
		$CG{$i+0} = 1 if $ref1 eq "C" and $ref2 eq "G";
		$CG{$i+1} = 1 if $ref1 eq "C" and $ref2 eq "G";
		if ($ref1 eq "C") {
			push(@{$peak->{CG}{Conlyconv}}, $peak->{CG}{conv}[$i]);
			$CGind ++;
			$pos{C}{CGN}[$CGind] = $i;
			$pos{C}{NCG}[$i] = $CGind;
			$pos{C}{CGNseq}[$CGind] = "C";
			$pos{C}{NCGseq}[$i] = "C";
		}
		else {
			$pos{C}{NCGseq}[$i] = ".";
			$pos{C}{NCG}[$i] = $CGind;
		}
		if ($ref1 eq "C" and $ref2 ne "G") {
			push(@{$peak->{CH}{Conlyconv}}, $peak->{CH}{conv}[$i]);
			$CHind ++;
			$pos{C}{CHN}[$CHind] = $i;
			$pos{C}{NCH}[$i] = $CHind;
			$pos{C}{CHNseq}[$CHind] = "G";
			$pos{C}{NCHseq}[$i] = "G";
		}
		else {
			$pos{C}{NCHseq}[$i] = ".";
			$pos{C}{NCH}[$i] = $CHind;
		}
		if ($ref1 eq "G") {
			push(@{$peak->{GC}{Conlyconv}}, $peak->{GC}{conv}[$i]);
			$GCind ++;
			$pos{G}{GCN}[$GCind] = $i;
			$pos{G}{NGC}[$i] = $GCind;
			$pos{G}{GCNseq}[$GCind] = "G";
			$pos{G}{NGCseq}[$i] = "G";
		}
		else {
			$pos{G}{NGCseq}[$i] = ".";
			$pos{G}{NGC}[$i] = $GCind;
		}
		if ($ref0 ne "C") {
			push(@{$peak->{GH}{Conlyconv}}, $peak->{GH}{conv}[$i]);
			$GHind ++;
			$pos{G}{GHN}[$GHind] = $i;
			$pos{G}{NGH}[$i] = $GHind;
			$pos{G}{GHNseq}[$GHind] = "G";
			$pos{G}{NGHseq}[$i] = "G";
		}
		else {
			$pos{G}{NGHseq}[$i] = ".";
			$pos{G}{NGH}[$i] = $GHind;
		}
	}
	my @GCN = @{$pos{G}{GCN}};
	my @CGN = @{$pos{C}{CGN}};
	my @NGC = @{$pos{G}{NGC}};
	my @NCG = @{$pos{C}{NCG}};
	my @GCNseq = @{$pos{G}{GCNseq}};
	my @CGNseq = @{$pos{C}{CGNseq}};
	my @GNHseq = @{$pos{G}{GHNseq}};
	my @CHNseq = @{$pos{C}{CHNseq}};
	my @NGCseq = @{$pos{G}{NGCseq}};
	my @NCGseq = @{$pos{C}{NCGseq}};
	my @NGHseq = @{$pos{G}{NGHseq}};
	my @NCHseq = @{$pos{C}{NCHseq}};
	my @CpeakCH = @{$peak->{CH}{Conlyconv}};
	my @CpeakCG = @{$peak->{CG}{Conlyconv}};
	my @CpeakGH = @{$peak->{GH}{Conlyconv}};
	my @CpeakGC = @{$peak->{GC}{Conlyconv}};
	for (my $i = 0; $i < @{$peak->{CG}{Conlyconv}}; $i++) {
		$peak->{CG}{Conlyconv2}[$i] = $peak->{CG}{Conlyconv}[$i];
	}
	for (my $i = 0; $i < @{$peak->{GC}{Conlyconv}}; $i++) {
		$peak->{GC}{Conlyconv2}[$i] = $peak->{GC}{Conlyconv}[$i];
	}
	for (my $i = 0; $i < @{$peak->{CH}{Conlyconv}}; $i++) {
		$peak->{CH}{Conlyconv2}[$i] = $peak->{CH}{Conlyconv}[$i];
	}
	for (my $i = 0; $i < @{$peak->{GH}{Conlyconv}}; $i++) {
		$peak->{GH}{Conlyconv2}[$i] = $peak->{GH}{Conlyconv}[$i];
	}
	($peak->{CG}{Conlyconv}, $peak->{CG}{Conlycoor}, $peak->{CG}{convperc}) = get_peak_conly($peak->{CG}{Conlyconv}, $window, $threshold);
	($peak->{CG}{Conlyconv2}, $peak->{CG}{Conlycoor2}, $peak->{CG}{convperc2}) = get_peak_conly($peak->{CG}{Conlyconv2}, $window, 0.5);
	($peak->{GC}{Conlyconv}, $peak->{GC}{Conlycoor}, $peak->{GC}{convperc}) = get_peak_conly($peak->{GC}{Conlyconv}, $window, $threshold);
	($peak->{GC}{Conlyconv2}, $peak->{GC}{Conlycoor2}, $peak->{GC}{convperc2}) = get_peak_conly($peak->{GC}{Conlyconv2}, $window, 0.5);
	($peak->{CH}{Conlyconv}, $peak->{CH}{Conlycoor}, $peak->{CH}{convperc}) = get_peak_conly($peak->{CH}{Conlyconv}, $window, $threshold);
	($peak->{CH}{Conlyconv2}, $peak->{CH}{Conlycoor2}, $peak->{CH}{convperc2}) = get_peak_conly($peak->{CH}{Conlyconv2}, $window, 0.5);
	($peak->{GH}{Conlyconv}, $peak->{GH}{Conlycoor}, $peak->{GH}{convperc}) = get_peak_conly($peak->{GH}{Conlyconv}, $window, $threshold);
	($peak->{GH}{Conlyconv2}, $peak->{GH}{Conlycoor2}, $peak->{GH}{convperc2}) = get_peak_conly($peak->{GH}{Conlyconv2}, $window, 0.5);
	my $mergecoorprint;
	print "$name\n" if defined $mergecoorprint;
	print " CG\n" if defined $mergecoorprint;
	($peak->{CG}{Conlycoor}, $peak->{CG}{Conlyconvpercbed}) = mergecoor(
	 $peak->{CG}{Conlycoor}, $peak->{CG}{Conlycoor2}, 
	 $peak->{CG}{convperc}, $peak->{CG}{convperc2}, 
	 \@{$pos{C}{CGN}}, $window, $minLen,$mergecoorprint
	);
	print " GC\n" if defined $mergecoorprint;
	
	($peak->{GC}{Conlycoor}, $peak->{GC}{Conlyconvpercbed}) = mergecoor(
	 $peak->{GC}{Conlycoor}, $peak->{GC}{Conlycoor2}, 
	 $peak->{GC}{convperc}, $peak->{GC}{convperc2}, 
	 \@{$pos{G}{GCN}}, $window, $minLen,$mergecoorprint
	);
	
	print " CH\n" if defined $mergecoorprint;
	($peak->{CH}{Conlycoor}, $peak->{CH}{Conlyconvpercbed}) = mergecoor(
	 $peak->{CH}{Conlycoor}, $peak->{CH}{Conlycoor2}, 
	 $peak->{CH}{convperc}, $peak->{CH}{convperc2}, 
	 \@{$pos{C}{CHN}}, $window, $minLen,$mergecoorprint
	);
	
	print " GH\n" if defined $mergecoorprint;
	($peak->{GH}{Conlycoor}, $peak->{GH}{Conlyconvpercbed}) = mergecoor(
	 $peak->{GH}{Conlycoor}, $peak->{GH}{Conlycoor2}, 
	 $peak->{GH}{convperc}, $peak->{GH}{convperc2}, 
	 \@{$pos{G}{GHN}}, $window, $minLen,$mergecoorprint
	);

	my @peakCH = (".")x@NCHseq;
	my @peakCG = (".")x@NCGseq;
	my @peakGH = (".")x@NGHseq;
	my @peakGC = (".")x@NGCseq;
	$peak->{CH}{tot} = 0;
	$peak->{GH}{tot} = 0;
	$peak->{CG}{tot} = 0;
	$peak->{GC}{tot} = 0;

	foreach my $begC (sort {$a <=> $b} keys %{$peak->{CG}{Conlycoor}}) {
		my $endC = $peak->{CG}{Conlycoor}{$begC};
		my $beg = $pos{C}{CGN}[$begC+1];
		my $end = $pos{C}{CGN}[$endC+1];

		$peak->{CG}{peakcoor}{$beg} = $end;
		$peak->{CG}{convpercbed}{$beg} = $peak->{CG}{Conlyconvpercbed}{$begC};
		$peak->{CG}{tot} ++;

		for (my $i = $beg; $i < $end; $i++) {
			$peakCG[$i] = $peak->{CG}{tot};
			$peak->{CG}{conv}[$i] =~ tr/67/89/;
		}
	}

	foreach my $begC (sort {$a <=> $b} keys %{$peak->{CH}{Conlycoor}}) {
		my $endC = $peak->{CH}{Conlycoor}{$begC};
		my $beg = $pos{C}{CHN}[$begC+1];
		my $end = $pos{C}{CHN}[$endC+1];

		$peak->{CH}{peakcoor}{$beg} = $end;
		$peak->{CH}{convpercbed}{$beg} = $peak->{CH}{Conlyconvpercbed}{$begC};
		$peak->{CH}{tot} ++;
		for (my $i = $beg; $i < $end; $i++) {
			$peakCH[$i] = $peak->{CH}{tot};
			$peak->{CH}{conv}[$i] =~ tr/57/11/;
			$peak->{CH}{conv}[$i] =~ tr/6/8/;
		}
	}

	my @peakG = (".")x@NGCseq;
	foreach my $begG (sort {$a <=> $b} keys %{$peak->{GC}{Conlycoor}}) {
		my $endG = $peak->{GC}{Conlycoor}{$begG};
		my $beg = $pos{G}{GCN}[$begG+1];
		my $end = $pos{G}{GCN}[$endG+1];

		$peak->{GC}{peakcoor}{$beg} = $end;
		$peak->{GC}{convpercbed}{$beg} = $peak->{GC}{Conlyconvpercbed}{$begG};
		$peak->{GC}{tot} ++;
		for (my $i = $beg; $i < $end; $i++) {
			$peakGC[$i] = $peak->{GC}{tot};
			$peak->{GC}{conv}[$i] =~ tr/67/89/;
		}
	}
	foreach my $beg (sort keys %{$peak->{GC}{coor}}) {
		my $end = $peak->{GC}{coor}{$beg};
	}

	foreach my $begG (sort {$a <=> $b} keys %{$peak->{GH}{Conlycoor}}) {
		$peakGH ++;
		my $endG = $peak->{GH}{Conlycoor}{$begG};
		my $beg = $pos{G}{GHN}[$begG+1];
		my $end = $pos{G}{GHN}[$endG+1];

		$peak->{GH}{peakcoor}{$beg} = $end;
		$peak->{GH}{convpercbed}{$beg} = $peak->{GH}{Conlyconvpercbed}{$begG};
		$peak->{GH}{tot} ++;
		for (my $i = $beg; $i < $end; $i++) {
			$peakGH[$i] = $peak->{GH}{tot};
			$peak->{GH}{conv}[$i] =~ tr/57/11/;
			$peak->{GH}{conv}[$i] =~ tr/6/8/;
		}
	}
	@peakCH = @{$peak->{CH}{conv}};
	@peakCG = @{$peak->{CG}{conv}};
	@peakGH = @{$peak->{GH}{conv}};
	@peakGC = @{$peak->{GC}{conv}};
	
	return ($peak);
}
sub mergecoor {
	my ($coorhash1, $coorhash2, $convperc1, $convperc2, $posarr, $window, $minLen, $print) = @_;
	
	my @convperc1 = split("\t", $convperc1);
	my @convperc2 = split("\t", $convperc2);
	foreach my $beg1C (sort {$a <=> $b} keys %{$coorhash1}) {
		my $end1C = $coorhash1->{$beg1C};
		my $beg1 = $posarr->[$beg1C+1];
		my $end1 = $posarr->[$end1C+1];
		my $convperc1 = $convperc1[$beg1C];
		print "C:$LCY$beg1C\t$end1C$N\tN:$beg1\t$end1\t$convperc1[$beg1C]\tcoorhash1\t+\n" if defined $print;
	}
	foreach my $beg2C (sort {$a <=> $b} keys %{$coorhash2}) {
		my $end2C = $coorhash2->{$beg2C};
		my $beg2 = $posarr->[$beg2C+1];
		my $end2 = $posarr->[$end2C+1];
		my $convperc2 = $convperc2[$beg2C];
		print "C:$LCY$beg2C\t$end2C$N\tN:$beg2\t$end2\t$convperc2[$beg2C]\tcoorhash2\t+\n" if defined $print;
	}
	foreach my $beg2C (sort {$a <=> $b} keys %{$coorhash2}) {
		my $end2C = $coorhash2->{$beg2C};
		my $beg2 = $posarr->[$beg2C+1];
		my $end2 = $posarr->[$end2C+1];
		my $convperc2 = $convperc2[$beg2C];
		if (not defined $coorhash1->{$beg2C}) {
			$coorhash1->{$beg2C} = $end2C;
			$convperc1[$beg2C] = $convperc2[$beg2C];
			print "C:$LCY$beg2C\t$end2C$N\tN:$beg2\t$end2\t$convperc1[$beg2C]\tcoorhash2\t+ ${LCY}indiv$N\n" if defined $print;
		}
		else {
			$coorhash1->{$beg2C} = $end2C if $coorhash1->{$beg2C} < $end2C;
			$convperc1[$beg2C] = $convperc2[$beg2C] if $convperc1[$beg2C] < $convperc2[$beg2C];
			print "C:$LCY$beg2C\t$end2C$N\tN:$beg2\t$end2\t$convperc1[$beg2C]\tcoorhash2\t+ ${LGN}combined$N\n" if defined $print;
		}
	}
	my ($coorhash3, $convpercarr3) = mergepeakcoor($coorhash1, $convperc1);
	
	my $coorhash4;
	my $convhash4;
	foreach my $beg3C (sort {$a <=> $b} keys %{$coorhash3}) {
		my $end3C = $coorhash3->{$beg3C};
		my $beg3 = $posarr->[$beg3C+1];
		my $end3 = $posarr->[$end3C+1];
		my $convperc3 = $convpercarr3->[$beg3C];
		print "UDNEFINED beg3 of beg3C-end3C=$LCY$beg3C-$end3C$N\n" if not defined $beg3;
		if ($end3 - $beg3 >= $minLen or ($end3 - $beg3 >= $window and $convperc3 >= 0.5)) {
			$coorhash4->{$beg3C} = $end3C;
			$convhash4->{$beg3C} = $convpercarr3->[$beg3C];
			print "C:$LCY$beg3C\t$end3C$N\tN:$beg3\t$end3\t$convperc3\tCOORHASH3 ($LGN used$N)\n" if defined $print;
		}
		else {
			print "C:$LCY$beg3C\t$end3C$N\tN:$beg3\t$end3\t$convperc3\tCOORHASH3 (not used)\n" if defined $print;
		}
	}
	print "\n" if defined $print;
	return($coorhash4, $convhash4);
}
sub get_peak_conly {
	my ($convarr, $window, $threshold) = @_;
	my @conv = @{$convarr};
	my @peak;
	my %peakcoor;
	my $is_peak = 0;
	my %conv;
	my @convperc;
	my $convwindow = @conv - $window;
	for (my $i = 0; $i < @conv-$window; $i++) {
		my $chunk = join("", @conv[$i..($i+$window-1)]);
		my $conv = $chunk =~ tr/67/67/;
		my $perc = int(1000*$conv / $window+0.5)/1000;
		$convperc[$i] = $perc;
		$chunk =~ tr/67/89/ if ($perc >= $threshold);
		my ($beg, $end, $beg50, $end50);
		my $j0 = $i;
		my $j1 = $i + $window;
		for (my $j = $j0; $j < $j1; $j++) { 
			my $conv2 = $conv[$j];
			if ($perc >= $threshold and $conv2 =~ /[67]/) {
				$beg = $j if not defined $beg;
				$end = $j;
			}
			$conv2 =~ tr/67/89/ if $perc >= $threshold;
			$peak[$j] = $conv2 if not defined $peak[$j];
			if ($peak[$j] !~ /(89)/) {
				$peak[$j] = $conv2;
			}
		}
		if (defined $beg) {
			if (not defined $peakcoor{$beg}) {
				$peakcoor{$beg} = $end;
			}
			elsif ($peakcoor{$beg} < $end) {
				$peakcoor{$beg} = $end;
			}
			else {
			}
		}
	}
	my $convperc = join("\t", @convperc);
	my ($peakcoor2hash, $junk) = mergepeakcoor(\%peakcoor, \@convperc);
	my %peakcoor2;
	%peakcoor2 = %{$peakcoor2hash} if defined $peakcoor2hash;
	foreach my $beg (sort keys %peakcoor) {
		my $end = $peakcoor{$beg};
	}
	foreach my $beg (sort keys %peakcoor2) {
		my $end = $peakcoor2{$beg};
	}
	return(\@peak, \%peakcoor2, $convperc);
}


sub mergepeakcoor {
	my ($peakcoorhash, $convpercarr, $print) = @_;
	my %peakcoor = %{$peakcoorhash};
	my ($lastbeg, $lastend);
	my %peakcoor2;
	my @convperc = $convpercarr =~ /ARRAY/ ? @{$convpercarr} : split("\t", $convpercarr);
	my @convperc2;
	foreach my $beg (sort {$a <=> $b} keys %peakcoor) {
		my $end = $peakcoor{$beg};
		if (not defined $lastbeg) {
			$lastbeg = $beg;
			$lastend = $end;

			if (not defined $convperc2[$lastbeg]) {
				$convperc2[$lastbeg] = $convperc[$lastbeg];
				$convperc2[$beg] = $convperc[$beg];
			}
			if ($convperc2[$lastbeg] < $convperc[$lastbeg]) {
				$convperc2[$lastbeg] = $convperc[$lastbeg];
				$convperc2[$beg] = $convperc[$beg];
			}
			else {
				$convperc[$lastbeg] = $convperc2[$lastbeg];
				$convperc[$beg] = $convperc2[$beg];
			}
			print "1  $LGN$beg-$end$N: not defined lastbeg so lastbeg is $LGN$lastbeg-$lastend$N $convperc2[$lastbeg]\n" if defined $print;
		}
		elsif (($beg <= $lastend and $beg >= $lastbeg) or ($end <= $lastend and $end >= $lastbeg)) {
			print "2a $LGN$beg-$end$N: end $end >= $lastend and $lastbeg <= beg $beg <= $lastend, so lastend = end ($end)" if defined $print;
			$lastend = $end if $end >= $lastend;
			$lastbeg = $beg if $beg <= $lastbeg;
			if ($convperc2[$lastbeg] < $convperc[$beg]) {
				$convperc2[$lastbeg] = $convperc[$beg];
				$convperc2[$beg] = $convperc[$beg];
			}
			if ($convperc[$beg] < $convperc2[$lastbeg]) {
				$convperc[$beg] = $convperc2[$lastbeg];
				$convperc[$beg] = $convperc2[$lastbeg];
			}
			print " so $LGN$lastbeg-$lastend$N ($convperc2[$lastbeg])\n" if defined $print;
		}
		else {
			print "2b $LGN$beg-$end$N: NOT INTERSECT)\n" if defined $print;
			$peakcoor2{$lastbeg} = $lastend;
			print "3  $LGN$beg-$end$N: final coor $YW$lastbeg-$lastend$N $convperc2[$lastbeg]\n" if defined $print;
			$lastbeg = $beg;
			$lastend = $end;
			if (defined $convperc2[$beg]) {
				$convperc2[$lastbeg] = $convperc[$beg] if $convperc2[$lastbeg] < $convperc[$beg];
			}
			else {
				$convperc2[$lastbeg] = $convperc[$beg];
			}
			print "4  $LGN$beg-$end$N: new coor $YW$lastbeg-$lastend$N \n" if defined $print;
		}
	}
	if (defined $lastbeg) {
		print "5  endline and defined lastbeg so $LGN$lastbeg-$lastend$N\n" if defined $print;
		$peakcoor2{$lastbeg} = $lastend;
	}
	return(\%peakcoor2, \@convperc2);
}

sub det_peak2 { 
   my ($peak, $nuc, $type, $pos) = @_;
   $pos = "NA" if not defined $pos;
   return $peak if $peak ne '!';
   my @N =  $type eq 'CH' ? ('0', '0'   , 'cd', 'CDMNPQ', '12' , 'B\-' ) :
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
	for (my $i = 0; $i < $max; $i++) {
      my ($jmin, $jmax) = $i == 0 ? (0, $window) : ($i, $i+1);
      for (my $j = $jmin; $j < $jmax; $j++) {
         if ($j < $totalseqpos - $window) {
            if ($i == 0) {
               my $nucpos1 = $nuc->[$SEQ->{$gene}{loc}{pos1}{$j}];
               $con = det_CG_con($con, $nucpos1, 1);
               LOG($outLog, "i=$i, j=$j: nuc.i=$LCY$nucpos1$N, nuc.j(nucpos1)=$nucpos1: file.filtered.gz value isn't one of the listed at $LGN" . __LINE__ . "$N\n") if $nucpos1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
               $peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CH");
               $peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CG");#, $SEQ->{$gene}{loc}{pos1}{$j});
            }
            else {
               my ($nucpos0, $nucpos2) = ($nuc->[$SEQ->{$gene}{loc}{pos1}{$j-1}], $nuc->[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}]);
               $con = det_CG_con($con, $nucpos0, -1, $nucpos2, 1);
               $peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CH");
               $peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CG",$SEQ->{$gene}{loc}{pos1}{$j+$window-1});
            }
         }
         if ($j < $totalseqneg - $window) {
            if ($i == 0) {
               my $nucneg1 = $nuc->[$SEQ->{$gene}{loc}{neg1}{$j}];
               $con = det_CG_con($con, $nucneg1, 1);
               LOG($outLog, "i=$i, j=$j: nuc.i=$LCY$nucneg1$N, nuc.j(nucneg1)=$nucneg1: file.filtered.gz value isn't one of the listed at $LGN" . __LINE__ . "$N\n") if $nucneg1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
               $peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j}] = det_peak($peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j}], $nucneg1, "GH");
               $peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j}] = det_peak($peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j}], $nucneg1, "GC");
            }
            else {
               my ($nucneg0, $nucneg2) = ($nuc->[$SEQ->{$gene}{loc}{neg1}{$j-1}], $nuc->[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}]);
               $con = det_CG_con($con, $nucneg0, -1, $nucneg2, 1);
               $peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GH}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}], $nucneg2, "GH");
               $peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GC}[$SEQ->{$gene}{loc}{neg1}{$j+$window-1}], $nucneg2, "GC");
            }
         }
      }
      my ($CG, $CH, $GC, $GH) = (0,0,0,0);
      ($peak, $CG, $CH, $GC, $GH) = isPeak($i, $con, $peak, $window, $threshold, $SEQ->{$gene}, $totalseqpos, $totalseqneg, $check);
      $peakCH += $CH;
      $peakCG += $CG;
      $peakGC += $GC;
      $peakGH += $GH;
   }
   return ($peak, $peakCH, $peakCG, $peakGH, $peakGC);
}

sub det_peak { 
   my ($peak, $nuc, $type, $pos) = @_;
   $pos = "NA" if not defined $pos;

   return $peak if $peak ne '!';
   my @N =  $type eq 'CH' ? ('0', '0'   , 'cd', 'CDMNPQ', '12' , 'B\-' ) :
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

sub det_CG_con { 
   my ($con, @cons) = @_;
   for (my $i = 0; $i < @cons; $i+= 2) {
      my ($type, $val) = ($cons[$i], $cons[$i+1]);
      ($con->{CG}{con} += $val and $con->{CG}{tot} += $val and next) if $type eq "e";
      ($con->{GC}{con} += $val and $con->{GC}{tot} += $val and next) if $type eq "i";
      ($con->{CH}{con} += $val and $con->{CH}{tot} += $val and next) if $type eq "c" or $type eq "d";
      ($con->{GH}{con} += $val and $con->{GH}{tot} += $val and next) if $type eq "g" or $type eq "h";
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

   my $Ctotmin = 5;#$window;
   my $peakCG = ($CGtot >= $Ctotmin and ($CGcon / $CGtot) > $threshold) ? 1 : 0;
   my $peakCH = ($CHtot >= $Ctotmin and ($CHcon / $CHtot) > $threshold) ? 1 : 0;
   my $peakGC = ($GCtot >= $Ctotmin and ($GCcon / $GCtot) > $threshold) ? 1 : 0;
   my $peakGH = ($GHtot >= $Ctotmin and ($GHcon / $GHtot) > $threshold) ? 1 : 0;
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
            $temp =~ s/8([01345679A-Z\! ]+)$/2$1/i;
         }
         @{$peak->{CH}}[$min..$max] = split("", $temp);
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
            $temp =~ s/8([012345679A-Z\! \.]+)$/2$1/i;
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
   return ($peak, $peakCG, $peakCH, $peakGC, $peakGH);
}
sub mydefined {
   my @arr = @_;
   for (my $i = 0; $i < @arr; $i++) {
      (print "$i: mydefined failed\n" and return 1) if not defined $arr[$i];
   }
   return 0;
}
