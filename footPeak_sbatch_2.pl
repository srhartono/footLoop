#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);

use myFootLib;
use FAlite;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);

my ($thisfileName) = getFilename($0, "fullname");
my ($thismd5) = getMD5_simple($0);

my ($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $version_small) = @ARGV;
die "\nUsage: $YW$0$N ${LCY}indexFile$N ${LGN}seqFile$N ${LCY}i$N ${LGN}origFile$N ${LCY}outDir$N ${LGN}window$N ${LCY}threshold$N ${LGN}totalorigFile$N ${LCY}label$N ${LGN}genewant$N ${LCY}minDis$N ${LGN}minLen$N ${LCY}version_small$N\n\n" unless @ARGV == 13;
my ($origFolder, $origFilename) = getFilename($origFile,"folderfull");
undef $genewant if $genewant eq -1;
my $outLogFile = "$outDir/.footPeak_sbatch/$origFilename.footPeak_sbatch.log.txt";
open (my $outLog, ">", $outLogFile) or die "Failed to write to $LCY$outLogFile$N: $!\n";
my $SEQ = parse_indexFile_and_seqFile($indexFile, $seqFile, $outLog);

run_footPeak_old($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog);
main_2($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog);

close $outLog;
print "\n\n$outLogFile\n";
sub run_footPeak_old {
	#my $Q = new Thread::Queue;
	my ($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog) = @_;
	my $out;
	my $peakFile = $origFile;
	LOG($outLog, date() . "$LGN$origFileInd$N/$LGN$totalorigFile$N:$LPR footPeak.pl processing$N $YW$peakFile$N\n");

	if (defined $genewant and $peakFile !~ /$genewant/i) {
		LOG($outLog, date() . "-> Skipped $LCY$peakFile$N as it doesn't contain $LGN-G $genewant$N\n","NA");
		return;
	}
	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
	$peakFilename =~ s/.filtered.gz$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);

	
	DIELOG($outLog, "gene=$gene seq->gene{seq} isn't defined!\n") if not defined $SEQ->{$gene} or (defined $SEQ->{$gene} and not defined $SEQ->{$gene}{seq});

	my ($totalPeak, $linecount, $peakcount, $total, %data, %end, %bad, %final, %group) = (0,0,0, scalar(@{$SEQ->{$gene}{seq}}));
	#LOG($outLog, date() . "$YW$origFileInd/" . scalar(@origFile) . "$N DEBUG FILENAME=$peakFilename, GENE=$gene\n") if $origFileInd % 10 == 0;
	DIELOG($outLog, "Died cannot get gene from peakfile $peakFile\n") unless defined $gene;
	DIELOG($outLog, "Cannot find sequence of gene=$LCY$gene$N in $seqFile!\n") if not defined $SEQ->{$gene};
	my %pk = ("CH" => 0, "CG" => 0, "GH" => 0, "GC" => 0);

	#next if $size eq 0;

	foreach my $type (sort keys %pk) {
		my $outpeak = "PEAK$type";
		my $outnopk = "NOPK$type";
		open ($out->{$outpeak}, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$type.PEAK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_CG.PEAK: $!\n");
		open ($out->{$outnopk}, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$type.NOPK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_CG.NOPK: $!\n");
	}
	my ($peakFilesize) = -s $peakFile; #zgrep -c ^ $peakFile`; chomp($peakFilesize); # =~ /^(\d+)($|\s+$)/;
	LOG($outLog, "$LGN$origFileInd$N. size=$LGN$peakFilesize$N Gene = $LRD$gene$N total length = $LPR$total$N\n","NA") if $peakFilesize <= 20;
	#LOG($outLog, "$LGN$origFileInd$N. size=$LGN$peakFilesize$N Gene = $LRD$gene$N total length = $LPR$total$N\n") if $peakFilesize > 20;
	#next if $peakFilesize <= 20;

	my ($peakFilelinecount) = `zcat $peakFile|wc -l` =~ /^(\d+)$/;
	LOG($outLog, "$LGN$origFileInd$N/$LGN" . $totalorigFile .  "$N $LGN$peakFilename$N $LGN$peakFilelinecount$N size=$LGN$peakFilesize$N Gene = $LRD$gene$N total length = $LPR$total$N totalread=$LCY$peakFilelinecount$N\n") if $peakFilelinecount >= 5;
	#next if $peakFilelinecount < 5;
	#my ($l0, $t0) = (0,Benchmark->new());
	foreach my $type (sort keys %pk) {
		my $outpeak = "PEAK$type";
		my $outnopk = "NOPK$type";
		print {$out->{$outpeak}} "$peakFile\tPEAK\t$gene\t$type\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
		print {$out->{$outnopk}} "$peakFile\tNOPK\t$gene\t$type\t$strand\t" . join("\t", @{$SEQ->{$gene}{seq}}) . "\n";
	}


	my $in1;
	if ($peakFile =~ /.gz$/) {
		open ($in1, "zcat $peakFile|") or die "Cannot read from $peakFile: $!\n";
	}
	else {
		open ($in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
	}

	while (my $line = <$in1>) {
		chomp($line);
		LOG($outLog, date() . "${LPR}footPeak.pl$N: Done $LGN$linecount$N\n") if $linecount % 100 == 0;
		$linecount ++;
		next if $line =~ /^#/;
		my ($name, $type, $val) = split("\t", $line);
		$name = "footPeak.pl.name.UNDEF" if not defined $name or $name eq "" or $name =~ /^[ \t\n\s]*$/;
		$type = "footPeak.pl.type.UNDEF" if not defined $type or $type eq "" or $type =~ /^[ \t\n\s]*$/;
		$val = "footPeak.pl.val.UNDEF" if not defined $val or $val eq ""     or $val  =~ /^[ \t\n\s]*$/;
		if ($name eq "footPeak.pl.name.UNDEF" or $type eq "footPeak.pl.type.UNDEF" or $val eq "footPeak.pl.val.UNDEF") {
			LOG($outLog, "\n\n" . date() . "file=$YW$peakFilename$N line number $LGN$linecount$LRD ERROR$N: Skipped line error trying to parse$LCY tab delim$N :LGN\$name=$name $LPR\$type=$type $LCY\$val=$val$N, line:\n\n$line\n\n","NA");
			next;
		}
		$val = [split("", $val)];
		$name = "$gene.$name";
		my ($count, $beg, $end) = (0,0,0);
		my $check = 0;

		if (defined $genewant and $name !~ /$genewant/i) {$check = 1; next;} #genewant

		my ($pkcurr, $outcurr) = dothis($gene, $strand, $val, $name, $window, $threshold, $check);
		foreach my $type (sort keys %{$pkcurr}) {
			$pk{$type} += $pkcurr->{$type};
		}
		foreach my $type (sort keys %{$outcurr}) {
			print {$out->{$type}} $outcurr->{$type};
		}
		#last if $linecount >= 100;
	}

	my $to;
	foreach my $key (sort keys %pk) {
		next if $key =~ /NO$/;
		my $peak = defined $pk{$key} ? $pk{$key} : 0;
		$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
		LOG($outLog, "$key=$LGN$peak$N, ","NA");
	}
	LOG($outLog, "\n\n" . date() . "\n\n","NA");
	#foreach my $type (sort keys %pk) {
	#	my $outpeak = "PEAK$type";
	#	my $outnopk = "NOPK$type";
	#	close ($out->{$outpeak});#, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$type.PEAK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_CG.PEAK: $!\n");
	#	close ($out->{$outnopk});#, ">", "$outDir/.CALL/$peakFilename\_$window\_$threshold\_$type.NOPK") or DIELOG($outLog, "Failed to write to $outDir/.CALL/$peakFilename\_$window\_$threshold\_CG.NOPK: $!\n");
	#}
}

sub main_2 {
	my ($indexFile, $seqFile, $origFileInd, $origFile, $outDir, $window, $threshold, $totalorigFile, $label, $genewant, $minDis, $minLen, $SEQ, $version_small, $outLog) = @_;
	my $out;
	my $peakFile = $origFile;
	LOG($outLog, date() . "$LGN$origFileInd$N/$LGN$totalorigFile$N:$LPR footPeak.pl processing$N $YW$peakFile$N\n");

	if (defined $genewant and $peakFile !~ /$genewant/i) {
		LOG($outLog, date() . "-> Skipped $LCY$peakFile$N as it doesn't contain $LGN-G $genewant$N\n","NA");
		return;
	}

	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
	$peakFilename =~ s/.filtered.gz$//;
	$peakFilename = "$label\_gene$peakFilename";
	my ($gene, $strand) = $peakFilename =~ /_gene(.+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	DIELOG($outLog, "gene=$gene seq->gene{seq} isn't defined!\n") if not defined $SEQ->{$gene} or not defined $SEQ->{$gene}{seq};


   my $peakFilez = "$outDir/.CALL/$peakFilename\_$window\_$threshold\_CG.PEAK";
	my $footPeakres = main(($peakFilez, $seqFile, $gene, $minDis, $outDir, $minLen, $SEQ, $origFileInd, $totalorigFile, $version_small, $outLog));
	DIELOG($outLog, "\n" . date() . "Failed footpeak\n") if defined $footPeakres and $footPeakres eq -1;
}

sub main {

	# From footPeak.pl: 
	my ($input1, $faFile, $mygene, $minDis, $resDir, $minLen, $SEQ, $mygenecount, $totalgenecount, $version_small, $outLog) = @_;

#	LOG($outLog, "
#-----------------
#$YW $0 $version_small $N
#---------
#","NA");
#

	my ($OUTDIRS, $PEAK, $TEMP, $RCONV, $CPG, $ALL);
	($OUTDIRS->{FOOTPEAK}, $PEAK, $TEMP, $RCONV, $CPG, $ALL) = makeOutDir($resDir . "/.FOOTPEAK/");

	my @foldershort = split("\/", $resDir);
	my $foldershort = pop(@foldershort);
	($input1) = getFullpath($input1);
	my ($folder, $fileName) = getFilename($input1, "folderfull");

	makedir("$resDir/.CALL") if not -d "$resDir/\.CALL";
	makedir("$resDir/PEAKS_GENOME") if not -d "$resDir/PEAKS_GENOME";
	makedir("$resDir/PEAKS_LOCAL") if not -d "$resDir/PEAKS_LOCAL";

	my @coor = split("\t", $SEQ->{$mygene}{coor});
	my ($chr0, $beg0, $end0, $name0, $val0, $strand0) = @coor;
	die "Undefined beg or end at coor=\n" . join("\n", @coor) . "\n" if not defined $beg0 or not defined $end0;
	my $geneStrand = $strand0 eq "+" ? "Pos" : $strand0 eq "-" ? "Neg" : $strand0 =~ /^(Pos|Neg|Unk)$/ ? $strand0 : DIELOG($outLog, "Failed to parse gene strand from strand=$strand0, coor=$SEQ->{$mygene}{coor}\n");

	#LOG($outLog, "$mygenecount/$totalgenecount $mygene\n");
	my (%pk, %Rscripts, %files);
	my $total = make_total_hash();

	my $label = "";
	if (-e "$resDir/.LABEL") {
	   ($label) = `cat $resDir/.LABEL`;
	   chomp($label);
	}
	else {
		DIELOG($outLog, "Failed to parse label from .LABEL in $resDir/.LABEL\n");
	}
	my $parseName = parseName($fileName);
	my ($label2, $gene, $readStrand, $window, $thres, $rconvType) = @{$parseName->{array}};
	my $isPeak = $fileName =~ /\.PEAK/ ? "PEAK" : "NOPK";
	LOG($outLog, "\n\n-----" . date() . "$YW WARNING$N Inconsistent label in filename $LCY$fileName$N." . date() . " Label from $resDir/.LABEL: $label\nBut from fileName: $label2\n-----\n\n") if $label ne $label2;
	$label = $label2;

	if (defined $mygene and uc($mygene) ne uc($gene)) {
		LOG($outLog, date() . "$LRD ERROR$N: footPeakAddon.pm filename=$LCY$fileName$N, mygene=$mygene, input genez=$gene are not the same!\n") and return -1;
	}
	else {
		$mygene = $gene;
	}
	
	my $bad = find_lots_of_C($faFile, $mygene, $outLog) if defined $faFile;

	LOG($outLog, "\n\n------------$YW\n3. Peak calling$N\n","NA");
	LOG($outLog, date() . "$input1; Undefined mygene=$mygene=, strand=$readStrand=, window=$window=, thres=$thres=, type=$rconvType=, isPeak=$isPeak=\n") and exit 1 if not defined $isPeak or not defined $window;
	LOG($outLog, date() . "\n\nFolder $YW$folder$N\n","NA");
	my ($peakPrint, $nopkPrint, $headerPrint, $totalPrint);
	$headerPrint = "type\tconv\ttotal_read_unique_with_peak\ttotal_read_unique_without_peak\n";
	my @types = qw(CH CG GH GC);
	for (my $h = 0; $h < 4; $h++) {
		my $rconvType = $types[$h];
		my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK";
		my $nopkFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK";
		$files{$peakFile} = 1;
		LOG($outLog, "\n  ${LGN}3.$h$N ($rconvType)\n  - peakFile=$LCY$peakFile$N\n  - nopkFile=$LPR$nopkFile$N\n","NA");
	
		my ($folder1, $peakfileName) = getFilename($peakFile, "folderfull");
		my ($folder2, $nopkfileName) = getFilename($nopkFile, "folderfull");

		my $data;
		my ($linecount, $totalpeak, $totalnopk, $totalline) = (0,0,0,0);
		if (-e $nopkFile) {
			($totalline) = `wc -l $nopkFile` =~ /^\s*(\d+)/;
			$linecount = 0;
			open (my $in1, "<", $nopkFile) or LOG($outLog, date() . "Cannot read from $nopkFile: $!\n") and exit 1;
			LOG($outLog, date . " -> Processing NOPK file ($LGN$totalline$N lines)\n","NA");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . " --> Done $totalnopk / $totalline\n","NA") if $totalnopk % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, $bad, $minDis, $minLen, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data->{peak}}, $val) if $totalPeak > 0;
				push(@{$data->{nopk}}, $val) if $totalPeak == 0;
				$totalnopk ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in1;
		}
		my $peakCount = defined $data->{peak} ? @{$data->{peak}} : 0;
		my $nopkCount = defined $data->{nopk} ? @{$data->{nopk}} : 0;
		my $currflag;
		my $flag = getFlag($nopkFile, $geneStrand, $readStrand, $rconvType);
		$flag =~ s/^NOPK//;
		$total->{$rconvType}{peak}  += $peakCount;
		$total->{$rconvType}{nopk}  += $nopkCount;
		$total->{$rconvType}{total} += $totalnopk;
		my $currtotalline = $totalline;
		$nopkPrint .= "from_nopkFile$flag\t$rconvType\t$peakCount\t$nopkCount\t$totalnopk\n";

		if (-e $peakFile) {
			($totalline) = `wc -l $peakFile` =~ /^\s*(\d+)/;
			$linecount = 0;
			open (my $in1, "<", $peakFile) or LOG($outLog, date() . "Cannot read from $peakFile: $!\n") and exit 1;
			LOG($outLog, date . " -> Processing PEAK file ($LGN$totalline$N lines)\n","NA");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . " --> Done $totalpeak / $totalline\n","NA") if $totalpeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, $bad, $minDis, $minLen, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data->{peak}}, $val) if $totalPeak > 0;
				push(@{$data->{nopk}}, $val) if $totalPeak == 0;
				$totalpeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
				#last if $linecount == 100;
			}
			close $in1;
		}
		$peakCount = defined $data->{peak} ? @{$data->{peak}} - $peakCount : 0;
		$nopkCount = defined $data->{nopk} ? @{$data->{nopk}} - $nopkCount : 0;


		$total->{$rconvType}{peak}  += $peakCount;
		$total->{$rconvType}{nopk}  += $nopkCount;
		$total->{$rconvType}{total} += $totalpeak;
		$currtotalline += $totalline;
		$flag = getFlag($peakFile, $geneStrand, $readStrand, $rconvType);
		$flag =~ s/^PEAK//;
		$peakPrint  .= "from_peakFile$flag\t$rconvType\t$peakCount\t$nopkCount\n";
		#$totalPrint .= "PEAK$flag\t$rconvType\t$total->{$rconvType}{peak}\t$total->{$rconvType}{nopk}\n";
		my $totalread = $total->{$rconvType}{peak} + $total->{$rconvType}{nopk};
#		$totalPrint .= "${LPR}PEAK$flag$N=$LGN$total->{$rconvType}{peak}$N/$LGN$linecount$N";
		$totalPrint .= " ${LPR}PEAK$flag$N=$LGN$total->{$rconvType}{peak}$N/$LGN$totalread$N";
		if (defined $data->{peak}) {
			die if @{$data->{peak}} != $total->{$rconvType}{peak};
			open (my $out1, ">", "$resDir/.CALL/$peakfileName.out") or LOG($outLog, date() . "Cannot write to $peakfileName.out: $!\n") and exit 1;
			foreach my $val (sort @{$data->{peak}}) {
				print $out1 "$val\n";
			}
			close $out1;
		}
		if (defined $data->{nopk}) {
			die if @{$data->{nopk}} != $total->{$rconvType}{nopk};
			open (my $out1, ">", "$resDir/.CALL/$nopkfileName.out") or LOG($outLog, date() . "Cannot write to $nopkfileName.out: $!\n") and exit 1;
			foreach my $val (sort @{$data->{nopk}}) {
				print $out1 "$val\n";
			}
			close $out1;
		}		
		LOG($outLog, "\n","NA");
	}
	
	#LOG($outLog, "\n\n-----------\n$LGN SUCCESS!$N\n\n","NA");
	LOG($outLog, "- $LRD$readStrand$N: $totalPrint\n");
	#LOG($outLog, "${YW}Summary of total peaks in gene=$LCY$mygene$N Strand $LRD$readStrand$N:$N\n" . myFootLib::prettyPrint("$headerPrint$totalPrint\n") . "\n\n");
	#LOG($outLog, "\n\n-----------\n\n");
	

	foreach my $file (sort keys %pk) {
		my ($pk_filename) = getFilename($file, 'full');
		open (my $outPEAKS, ">", "$resDir/PEAKS_GENOME/$pk_filename.genome.bed")   or LOG($outLog, "\tFailed to write into $resDir/PEAKS_GENOME/$pk_filename.genome.bed: $!\n")  and exit 1;
		open (my $outRPEAKS, ">", "$resDir/PEAKS_LOCAL/$pk_filename.local.bed") or LOG($outLog, "\tFailed to write into $resDir/PEAKS_LOCAL/$pk_filename.local.bed: $!\n") and exit 1;
		my $currtype = $file =~ /_CH/ ? "CH" : $file =~ /_CG/ ? "CG" : $file =~ /_GH/ ? "GH" : $file =~ /_GC/ ? "GC" : "UNK";
		LOG($outLog, "\tFailed to determine type (CH/CG/GH/GC) of file=$file.PEAKS: $!\n") if $currtype eq "UNK";
		foreach my $name (sort keys %{$pk{$file}}) {
			foreach my $peak (sort @{$pk{$file}{$name}}) {
				my ($beg, $end) = split("-", $peak);
				my $end1 = $end + $beg0;
				my $beg1 = $beg + $beg0;
				my ($junk, $readname) = $name =~ /^(\w+)\.(.+)$/;
				print $outPEAKS  "$chr0\t$beg1\t$end1\t$name\t0\t$strand0\t$file\n";
				print $outRPEAKS "$name\t$beg\t$end\n";
			}
		}
		close $outPEAKS;
		close $outRPEAKS;
	}
	
	makedir("$resDir/99_FOOTSTATS/") if not -d "$resDir/99_FOOTSTATS/";
	makedir("$resDir/99_FOOTSTATS/.PEAKSTATSTEMP") if not -d "$resDir/99_FOOTSTATS/.PEAKSTATSTEMP";
	open (my $outLGENE, ">", "$resDir/99_FOOTSTATS/.PEAKSTATSTEMP/.0_RESULTS\_$label\_gene$mygene\_$readStrand\_$window\_$thres.TXT");
	print $outLGENE "#folder\tpeak_file\tgene\tread_strand\tread_unique_total\tread_unique_with_peak_total\tread_unique_with_peak_perc\tpeakType\n";
	for (my $h = 0; $h < 4; $h++) {
		my $rconvType = $types[$h];
		my $totalPeak = $total->{$rconvType}{peak};
		my $totalNopk = $total->{$rconvType}{nopk};
		$total->{$rconvType}{peak} = $total->{$rconvType}{total} == 0 ? 0 : int(1000 * $total->{$rconvType}{peak} / $total->{$rconvType}{total}+0.5)/10;
		$total->{$rconvType}{nopk} = $total->{$rconvType}{total} == 0 ? 0 : int(1000 * $total->{$rconvType}{nopk} / $total->{$rconvType}{total}+0.5)/10;
		my @folder = split("/", $resDir);
		my $foldershort = $folder[@folder-1];
		   $foldershort = $folder[@folder-2] if not defined ($foldershort) or (defined $foldershort and $foldershort =~ /^[\s]*$/);
		my $peakFile    = "$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK";
		my $flag = getFlag($peakFile, $geneStrand, $readStrand, $rconvType);
		print $outLGENE "$foldershort\t$peakFile\t$mygene\t$rconvType\t$total->{$rconvType}{total}\t$totalPeak\t$total->{$rconvType}{peak}\t$flag\n";
	}
	close $outLGENE;

	my $stats = `cat $resDir/99_FOOTSTATS/.PEAKSTATSTEMP/.0_RESULTS\_$label\_gene$mygene\_$readStrand\_$window\_$thres.TXT`;
	my ($sampleName) = $folder =~ /^.+(PCB[\_\-]*\d+)/i;
	if ($folder =~ /debarcode/) {
		my ($temp) = $folder =~ /_ccs_(\w+)/;
		$sampleName = $sampleName . "_$temp";
	}
	if (not defined $sampleName) {
		$sampleName = $foldershort;
	}
	#LOG($outLog, "\n--------- Back to footPeak.pl --------\n\n\n","NA");
}

###############
# Subroutines #
###############


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
	#my $name_want = "PFC9.m64069_230413_005848/106890704/ccs";#AIRN_PFC66_FORWARD.16024";#CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/16024/ccs";
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
					LOG($outLog, date() . "\t\t$LGN YES$N peak=$beg-$end, bad=$begBad-$endBad\n","NA") if $isPeak eq "PEAK";# if $name eq "$name_want");
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
   my ($gene, $strand, $val, $name, $window, $threshold, $check) = @command;
   print "gene=$gene, strand=$strand, name=$name, Undefined seq!\n" and return 1 if not defined $SEQ;
   my ($peak, $peakCH, $peakCG, $peakGH, $peakGC) = getPeak(@command);
   my ($CH, $CG, $GH, $GC) = ("","","","");
   for (my $i = 0; $i < @{$SEQ->{$gene}{seq}}; $i++) {
      die "Undefined i=$LGN$i$N \$val->[\$i]\n" if not defined $val->[$i];
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
   return($pk, $out);
   #return 0;
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
#  return($peak,0,0,0,0); #debug
   for (my $i = 0; $i < $max; $i++) {
      my ($jmin, $jmax) = $i == 0 ? (0, $window) : ($i, $i+1);
      for (my $j = $jmin; $j < $jmax; $j++) {
         if ($j < $totalseqpos - $window) {
            if ($i == 0) {
               my $nucpos1 = $nuc->[$SEQ->{$gene}{loc}{pos1}{$j}];
               $con = det_CG_con($con, $nucpos1, 1);
               LOG($outLog, "i=$i, j=$j: nuc.i=$LCY$nucpos1$N, nuc.j(nucpos1)=$nucpos1: file.filtered.gz value isn't one of the listed at $LGN" . __LINE__ . "$N\n") if $nucpos1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
               $peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CH");
#              $peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CG");
               $peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j}], $nucpos1, "CG", $SEQ->{$gene}{loc}{pos1}{$j});
            }
            else {
               my ($nucpos0, $nucpos2) = ($nuc->[$SEQ->{$gene}{loc}{pos1}{$j-1}], $nuc->[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}]);
               $con = det_CG_con($con, $nucpos0, -1, $nucpos2, 1);
               $peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CH}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CH");
               $peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CG}[$SEQ->{$gene}{loc}{pos1}{$j+$window-1}], $nucpos2, "CG",$SEQ->{$gene}{loc}{pos1}{$j+$window-1});
            }
#           $chunkC .= $nucpos1;
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
#           $chunkG .= $nucneg1;
         }
#=cut
      }
      if ($i == 0 and defined $check) {
#        print "i=$i, " . join("", @{$peak->{CH}}[0..$window-1]) . "\n";
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
#  return 1; #debug
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

sub det_CG_con { # time= 5 line/s (57 l/s -> 52 l/s)
   my ($con, @cons) = @_;
#  return $con; #debug
   for (my $i = 0; $i < @cons; $i+= 2) {
      my ($type, $val) = ($cons[$i], $cons[$i+1]);
      ($con->{CG}{con} += $val and $con->{CG}{tot} += $val and next) if $type eq "e";
      ($con->{GC}{con} += $val and $con->{GC}{tot} += $val and next) if $type eq "i";
      ($con->{CH}{con} += $val and $con->{CH}{tot} += $val and next) if $type eq "c" or $type eq "d";
      ($con->{GH}{con} += $val and $con->{GH}{tot} += $val and next) if $type eq "g" or $type eq "h";
#     ($con->{CG}{tot} += $val and next) if $type eq "F" or $type eq "E" or $type eq "e" or $type eq "O" or $type eq "3" or $type eq "R";
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

   my $Ctotmin = 2;
   #instead of 5, becomes 2
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
            $temp =~ s/8([01345679A-Z\! ]+)$/2$1/i;# if $temp =~ /8([012345679A-Z\! ]+)$/i;
         }
         @{$peak->{CH}}[$min..$max] = split("", $temp);
         #die "$temp\n" if $peakCH eq 1;
   #     print "$i: $temp\n" if defined $check;
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
#     print "i=$i, min=$min, max=$max, CHcon=$LGN$CHcon$N, CHtot=$LPR$CHtot$N, " . join("", @{$peak->{CH}}[$min..$max]) . "\n" if defined $check and $check eq 1;# and $peakCH eq 1;
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
#  return ($peak, 0,0,0,0);
   return ($peak, $peakCG, $peakCH, $peakGC, $peakGH);
}
sub mydefined {
   my @arr = @_;
   for (my $i = 0; $i < @arr; $i++) {
      (print "$i: mydefined failed\n" and return 1) if not defined $arr[$i];
   }
   return 0;
}
