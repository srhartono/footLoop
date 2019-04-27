package footPeakAddon;

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_x $opt_R $opt_c $opt_t);

use myFootLib;
use FAlite;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);

my ($thisfileName) = getFilename($0, "fullname");
my ($thismd5) = getMD5_simple($0);

sub main {
	# From footPeak.pl: 
	my ($input1, $faFile, $mygene, $minDis, $resDir, $minLen, $SEQ, $version_small, $outLog) = @_;

	LOG($outLog, "
-----------------
$YW $0 $version_small $N
---------
");


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

	LOG($outLog, "\n\n------------$YW\n3. Peak calling$N\n\n\n");
	LOG($outLog, date() . "$input1; Undefined mygene=$mygene=, strand=$readStrand=, window=$window=, thres=$thres=, type=$rconvType=, isPeak=$isPeak=\n") and exit 1 if not defined $isPeak or not defined $window;
	LOG($outLog, date() . "\n\nFolder $YW$folder$N\n");
	my ($peakPrint, $nopkPrint, $headerPrint, $totalPrint);
	$headerPrint = "type\tconv\ttotal_read_unique_with_peak\ttotal_read_unique_without_peak\n";
	my @types = qw(CH CG GH GC);
	for (my $h = 0; $h < 4; $h++) {
		my $rconvType = $types[$h];
		my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK";
		my $nopkFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK";
		$files{$peakFile} = 1;
		LOG($outLog, "\n  ${LGN}3.$h$N ($rconvType)\n  - peakFile=$LCY$peakFile$N\n  - nopkFile=$LPR$nopkFile$N\n");
	
		my ($folder1, $peakfileName) = getFilename($peakFile, "folderfull");
		my ($folder2, $nopkfileName) = getFilename($nopkFile, "folderfull");

		my $data;
		my ($linecount, $totalpeak, $totalnopk, $totalline) = (0,0,0,0);
		if (-e $nopkFile) {
			($totalline) = `wc -l $nopkFile` =~ /^\s*(\d+)/;
			$linecount = 0;
			open (my $in1, "<", $nopkFile) or LOG($outLog, date() . "Cannot read from $nopkFile: $!\n") and exit 1;
			LOG($outLog, date . " -> Processing NOPK file ($LGN$totalline$N lines)\n");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . " --> Done $totalnopk / $totalline\n") if $totalnopk % 500 == 0;
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
			LOG($outLog, date . " -> Processing PEAK file ($LGN$totalline$N lines)\n");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . " --> Done $totalpeak / $totalline\n") if $totalpeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, $bad, $minDis, $minLen, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data->{peak}}, $val) if $totalPeak > 0;
				push(@{$data->{nopk}}, $val) if $totalPeak == 0;
				$totalpeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
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
		$totalPrint .= "PEAK$flag\t$rconvType\t$total->{$rconvType}{peak}\t$total->{$rconvType}{nopk}\n";

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
		LOG($outLog, "\n");
	}
	
	LOG($outLog, "\n\n-----------\n$LGN SUCCESS!$N\n\n");
	LOG($outLog, "${YW}Summary of total peaks in gene=$LCY$mygene$N Strand $LRD$readStrand$N:$N\n" . myFootLib::prettyPrint("$headerPrint$totalPrint\n") . "\n\n");
	LOG($outLog, "\n\n-----------\n\n");
	

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
	LOG($outLog, "\n--------- Back to footPeak.pl --------\n\n\n");
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
	my $name_want = "AIRN_PFC66_FORWARD.16024";#CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/16024/ccs";
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
	LOG($outLog, "\n------------\n${YW}2. Getting polyC (or G) regions from sequence file $LCY$seqFile$N\n\n\n");
	open(my $SEQIN, "<", $seqFile) or LOG($outLog, date() . "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!") and exit 1;
	my $fasta = new FAlite($SEQIN);
	my %lotsOfC;

	while (my $entry = $fasta->nextEntry()) {
	   my $gene = uc($entry->def); $gene =~ s/^>//;
	   my $seqz = uc($entry->seq);
		next if $gene ne $mygene;
	   LOG($outLog, date . "gene=$LGN$gene$N\n");
	   
	
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
	      LOG($outLog, date() . "$gene: $beg_C to $end_C ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n","NA");
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
	      LOG($outLog, date() . "$gene: $beg_G to $end_G ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n","NA");
	      $lotsOfC{$gene} .= "G;$beg_G;$end_G,";
	   }
		my $bad;
		if (defined $lotsOfC{$gene}) {
			LOG($outLog, date() . "Region of Bad C's:\n");
			open (my $outBadC, ">", "$seqFile.badc.bed") or DIELOG($outLog, "Failed to write to $seqFile.badc.bed: $!\n");
			$lotsOfC{$gene} =~ s/,$//;
			my @lotsOfC = split(",", $lotsOfC{$gene});
			foreach my $coor (@lotsOfC) {
				my ($nuc, $beg, $end) = split(";", $coor);
				my $strands = $nuc eq "C" ? 0 : $nuc eq "G" ? 16 : 255;
				my $strandsymbol = $strands eq 0 ? "+" : "-";
				$bad->{$strands}{$beg} = $end-1;
				LOG($outLog, date() . "$strands: coor=$coor, nuc=$nuc, beg=$beg, end=$end, beg-end=$beg-$bad->{$strands}{$beg}\n");
				print $outBadC "$gene\t$beg\t$end\t$coor\t" . ($end-$beg) . "\t$strandsymbol\n";
			}
			LOG($outLog, "\n${LGN}SUCCESS!$N PolyC (or G) bed file:\n$LCY$seqFile.badc.bed$N\n\n");
			return $bad;
		}
	}
	return;
}

1;
