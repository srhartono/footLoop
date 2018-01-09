package footPeakAddon;

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_x $opt_R $opt_c $opt_t);

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";

sub main 
{
	my ($input1, $faFile, $mygene) = @_;
	die "\nusage: $YW$0$N [-c to use cpg] $CY<CALM3_Pos_20_0.65_CG.PEAK>$N $CY<location with lots of C>$N\n\n" unless @_ == 3;
	my %pk;
	die "Input cannot be directry!\n" if -d $input1;
	($input1) = getFullpath($input1);
	#inputs END1
	my ($folder, $fileName) = getFilename($input1, "folderfull");
	my ($genez, $strand, $window, $thres, $type, $isPeak) = $fileName =~ /^(\w+)_(Pos|Neg)_(\d+)_(\d+\.?\d*)_(\w+)\.(PEAK|NOPK)$/;
	$mygene = $genez if not defined $mygene;
	
	open (my $outLog, ">", "$folder/footLoop_addition_logFile.txt") or die;
	my $lotsOfC = find_lots_of_C($faFile, $mygene, $outLog) if defined $faFile;
	my %bad;
	if (defined $lotsOfC) {
		my @lotsOfC = split(",", $lotsOfC);
		foreach my $coor (@lotsOfC) {
			my ($nuc, $beg, $end) = split(";", $coor);
			my $strand = $nuc eq "C" ? 0 : $nuc eq "G" ? 16 : 255;
			$bad{$strand}{$beg} = $end-1;
		}
	}


	die "$input1; Undefined mygene=$mygene=, strand=$strand=, window=$window=, thres=$thres=, type=$type=, isPeak=$isPeak=\n" if not defined $isPeak or not defined $window;
	my %total; 
	$total{Pos}{peak} = 0; $total{Pos}{nopeak} = 0; $total{Pos}{total} = 0;
	$total{Neg}{peak} = 0; $total{Neg}{nopeak} = 0; $total{Neg}{total} = 0;
	$total{Unk}{peak} = 0; $total{Unk}{nopeak} = 0; $total{Unk}{total} = 0;
	my $log2 = "";
	print STDERR "\n\nFolder $YW$folder$N: Processing files related to $LCY$input1$N\n";
	open (my $outLGENE, ">", "$folder/.0_RESULTS\_$mygene.TXT") if not defined $opt_x;
	open (my $outLEXTRA, ">", "$folder/.1_RESULTS_EXTRA\_$mygene.TXT") if not defined $opt_x;
	for (my $h = 0; $h < 4; $h++) {
		$type = $h % 4 == 0 ? 'CH' : $h % 4 == 1 ? 'CG' : $h % 4 == 2 ? 'GH' : 'GC';
		my $peakFile   = "$folder/$mygene\_$strand\_$window\_$thres\_$type.PEAK";
		my $nopkFile   = "$folder/$mygene\_$strand\_$window\_$thres\_$type.NOPK";
		print STDERR "h=$LGN$h\t$YW$peakFile\t$LCY$nopkFile\n$N";
	
		my ($folder1, $peakfileName) = getFilename($peakFile, "folder");
		my ($folder2, $nopeakfileName) = getFilename($nopkFile, "folder");
		
		my %data; my $totalnopeak = 0;
		my $linecount = 0;
		my $totalline = 0;
		if (-e $nopkFile) {
			($totalline) = `wc -l $nopkFile` =~ /^(\d+)/;
			open (my $in2, "<", $nopkFile) or die "Cannot read from $nopkFile: $!\n";
			print STDERR "\tProcessing $LPR$nopkFile$N ($LGN$totalline$N lines)\n";
			while (my $line = <$in2>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1;
				print STDERR "\tDone $totalnopeak / $totalline\n" if $totalnopeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, \%bad);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data{peak}}, $val) if $totalPeak > 0;
				push(@{$data{nopeak}}, $val) if $totalPeak == 0;
				$totalnopeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in2;
		}
	
		my $peakCount = defined $data{peak} ? @{$data{peak}} : 0;
		my $nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} : 0;
		my $nopeakPrint ="$folder2\t$nopeakfileName\t$peakCount\t$nopeakCount\t$totalnopeak\t$totalline";
		$total{$type}{peak} += $peakCount;
		$total{$type}{nopeak} += $nopeakCount;
		$total{$type}{total} += $totalnopeak;
		
		my $totalpeak = 0;
		if (-e $peakFile) {
			open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
			($totalline) = `wc -l $peakFile` =~ /^(\d+)/;
			print STDERR "\tProcessing $LPR$peakFile$N ($LGN$totalline$N lines)\n";
			$linecount = 0;
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1;
				print STDERR "\tDone $totalpeak / $totalline\n" if $totalpeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, \%bad);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data{peak}}, $val) if $totalPeak > 0;
				push(@{$data{nopeak}}, $val) if $totalPeak == 0;
				$totalpeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in1;
		}
		$peakCount = defined $data{peak} ? @{$data{peak}} - $peakCount : 0;
		$nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} - $nopeakCount : 0;
		my $peakPrint ="$folder1\t$peakfileName\t$peakCount\t$nopeakCount\t$totalpeak\t$totalline";
		$total{$type}{peak} += $peakCount;
		$total{$type}{nopeak} += $nopeakCount;
		$total{$type}{total} += $totalpeak;
	
	
	
		open (my $out1, ">", "$folder1/$peakfileName.out") or die "Cannot write to $peakfileName.out: $!\n";
		open (my $out2, ">", "$folder1/$nopeakfileName.out") or die "Cannot write to $nopeakfileName.out: $!\n";
		if (defined $data{peak}) {
			foreach my $val (sort @{$data{peak}}) {
				print $out1 "$val\n";
			}
		}
		if (defined $data{nopeak}) {
			foreach my $val (sort @{$data{nopeak}}) {
				print $out2 "$val\n";
			}
		}
		
		close $out1;
		close $out2;
#### HERE ##
#	next;
	
		$log2 .= "\#Folder\tFile\tPeak\tNoPeak\tTotalRead\tTotalLineInFile\n" if $h == 0 and not defined $opt_x;
		print STDERR "$peakPrint\n$nopeakPrint\n" if defined $opt_x;
		$log2 .= "$peakPrint\n$nopeakPrint\n" if not defined $opt_x;
	
		mkdir "$folder1/remove" if not -d "$folder1/remove/";
		
		my $peakFileBackup = "$folder1/remove/$peakfileName";
		my $peakFileTemp   = $peakFileBackup;
		my $count = 0;
		while (-e $peakFileTemp) {
			$peakFileTemp = $peakFileBackup . $count;
			$count ++;
		}
		if (-e $peakFile) {
	#		print STDERR "\tmv $peakFile $peakFileTemp\n";
	#		system("/bin/mv $peakFile $peakFileTemp") if not defined $opt_x;
	#		print STDERR "\tmv $folder1/$peakfileName.out $folder1/$peakfileName\n";
	#		system("mv $folder1/$peakfileName.out $folder1/$peakfileName") if not defined $opt_x;
		}
		
		my $nopkFileBackup = "$folder1/remove/$nopeakfileName";
		my $nopkFileTemp   = $nopkFileBackup;
		$count = 0;
		while (-e $nopkFileTemp) {
			$nopkFileTemp = $nopkFileBackup . $count;
			$count ++;
		}
		if (-e $nopkFile) {
	#		print STDERR "\t/bin/mv $nopkFile $nopkFileTemp\n";
	#		system("/bin/mv $nopkFile $nopkFileTemp") if not defined $opt_x;
	#		print STDERR "\tmv $folder1/$nopeakfileName.out $folder1/$nopeakfileName\n";
	#		system("mv $folder1/$nopeakfileName.out $folder1/$nopeakfileName") if not defined $opt_x;
		}
	}
	#		foreach my $peak (sort @{$peak{peak}}) {
	#			my ($beg, $end) = split("-", $peak);
	
	foreach my $file (sort keys %pk) {
		next if not defined $pk{$file};
		next if defined $pk{$file} and keys %{$pk{$file}} == 0;
		open (my $outR, ">", "$file.PEAKS") or die;
		foreach my $name (sort keys %{$pk{$file}}) {
			next if not defined $pk{$file}{$name};
			my @arr = @{$pk{$file}{$name}};
			foreach my $peakz (sort @{$pk{$file}{$name}}) {
				my ($beg, $end) = split("-", $peakz);
				my $name2 = $name; $name2 =~ s/^SEQ_//;
				print $outR "$name2\t$beg\t$end\n";
			}
		}
		close $outR;
	}
	
	if (not defined $opt_x) {
		my @types = qw(CH CG GH GC);
		foreach my $types (sort @types) {
			$total{$types}{total} = 0 if not defined $total{$types}{total};
			$total{$types}{peak} = 0 if not defined $total{$types}{peak};
			$total{$types}{nopeak} = 0 if not defined $total{$types}{nopeak};
			$total{$types}{peak} = $total{$types}{total} == 0 ? 0 : int(1000 * $total{$types}{peak} / $total{$types}{total}+0.5)/10;
			$total{$types}{nopeak} = $total{$types}{total} == 0 ? 0 : int(1000 * $total{$types}{nopeak} / $total{$types}{total}+0.5)/10;
			my @folder = split("/", $folder);
			my $foldershort = $folder[@folder-1];
			   $foldershort = $folder[@folder-2] if not defined ($foldershort) or (defined $foldershort and $foldershort =~ /^[\s]*$/);
			print $outLGENE "#folder\tfilename\tGene\tStrand\ttotal\tpeak.perc\n";
			print $outLGENE "$foldershort\t$fileName\t$mygene\t$types\t$total{$types}{total}\t$total{$types}{peak}\n";
			print $outLEXTRA "$log2";
		}
		close $outLGENE;
		close $outLEXTRA;
	}
	system("cat $folder/.0_RESULTS\_$mygene.TXT") if not defined $opt_x;
	print STDERR "\tcd $folder && run_Rscript.pl *MakeHeatmap.R\n";
	system("cd $folder && run_Rscript.pl *MakeHeatmap.R") if not defined $opt_x and defined $opt_R;
}
###############
# Subroutines #
###############

sub parse_peak {
	my $bad = $_[1];
	my %bad; %bad = %{$bad} if defined $bad;
	my ($name, $isPeak, $mygene, $type, $strand, @val) = split("\t", $_[0]);
	my $name_want = "CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/9477/ccs";
	shift(@val) if $val[0] eq "";
#	for (my $i = 0; $i < @val; $i++) {
#		if ($val[$i] !~ /^[456789]$/) {print "."} else {print "$val[$i]";}
#		print "\n" if $i != 0 and ($i+1) % 100 == 0;
#	}
#	print "\n";
	my $peaks;
	my %peak; $peak{curr} = 0; my $edge = 0; my $edge2 = 0; my $zero = 0; my $edge1 = 0;
	my $Length = @val; 
	my $print = "Total length = $Length\n";
	for (my $i = 0; $i < @val; $i++) {
		if ($i % 100 == 0 and $edge2 <= 1) {$print .= "\n$YW" . $i . "$N:\t";}
		my $val = $val[$i];
		$zero ++ if $val == 0;
		$edge = 1 if $zero > 20;
		$zero = 0 if $val != 0;
		if ($edge == 1 and $zero == 0) {
			$edge = $i;
			$edge1 = $i;
			$print .= "EDGE1 = $edge\n";
			$edge2  = 1;
		}
		if ($edge2 == 1 and $zero >= 20) {
			$edge2 = $i - $zero;
			$print .= "EDGE2 = $edge2\n";
		}
		if ($val =~ /[89]/) { # Peak Converted CpG or CH
			$peak{curr} = 1;
			$print .= "${LRD}$val$N";
			if (not defined $peak{beg}) {
				$peak{beg} = $i;
				$peak{end} = $i;
			}
			elsif (defined $peak{beg} and $i - $peak{end} >= 250) {
				push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
				$peak{beg} = $i;
				$peak{end} = $i;
			}
			elsif (defined $peak{beg} and $i - $peak{end} < 250) {
				$peak{end} = $i;
			}
		}
		else {
			$print .= "${LGN}$val$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val$N" if $val =~ /^[57]$/;
			$print .= "." if $val eq 1;
			$print .= "x" if $val eq 0 and ($edge2 <= 1 or ($edge2 > 1 and $i < $edge2));
			$print .= "EDGE2\n" if $val eq 0 and ($edge2 > 1 and $i == $edge2);
			if ($peak{curr} == 1 and $i - $peak{end} >= 250) {
				push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
				undef $peak{beg}; undef $peak{end};
				$peak{curr} = 0;
			}
		}
	}
	if ($peak{curr} == 1 and defined $peak{beg}) {
		push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
	}
	my (%nopeak, @peak);
	my %peak2;
	$print .= "\n";
	print "\nDoing $YW$name$N\n" if $name eq $name_want;#"SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
	if (defined $peak{peak}) {
		foreach my $peak (sort @{$peak{peak}}) {
			my ($beg, $end) = split("-", $peak);
			print "$name: $beg to $end\n" if $name eq 77011 or $name eq "$name_want";
			my $checkBad = 0;
			foreach my $begBad (sort keys %bad) {
				my $endBad = $bad{$begBad};
				print "\t$beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "$name_want";
				if ($beg >= $begBad and $beg <= $endBad and $end >= $begBad and $end <= $endBad) {
#					for (my $m = $beg; $m <= $end; $m++) {
#						$nopeak{$m} = 1;
#					}
					print "\t\t$LGN YES$N\n" if $name eq "$name_want";
					$checkBad = 1; last;
				}
			}
			if ($checkBad != 1) {
				foreach my $begBad (sort keys %bad) {
					my $endBad = $bad{$begBad};
					print "$name: $beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "$name_want";
					#if (($beg >= $begBad and $beg <= $endBad) or ($end >= $begBad and $end <= $endBad)) {
						my @valz = @val;
						my ($goodC, $badC) = (0,0);
						for (my $m = $beg; $m < $end; $m++) {
							if ($m >= $begBad and $m <= $endBad) {
								$badC ++ if $valz[$m] =~ /^(8|9)$/;
							}
							else {
								$goodC ++ if $valz[$m] =~ /^(8|9)$/;
							}
						}
						if ($goodC < 5 and $badC >= 9) {
#							print "\t$LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
							print "\t$YW$name$N $LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n";
							$checkBad = 1; last;
						}
#						else {
#							#print "\t$LGN OKAY!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
#						}
					print "\t\t$LGN YES$N\n" if $name eq "$name_want";
					#}
				}
			}
#			print "\t$name checkbad = $checkBad\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
			
			if ($checkBad == 1) {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n";
				print "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n" if $name eq "$name_want";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
			elsif ($end > 10 + $edge1 and $beg < $edge2 - 10) {
				print "something wrong\n";# if $name eq "$name_want";
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Used: $LGN$peak$N\n";
				push(@peak, "$beg-$end");
				push(@{$peak2{peak}}, $peak);
			}
			else {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
		}
	}
	my $totalpeak = scalar(@peak);
#	my @val2;
#	for (my $i = 0; $i < @val; $i++) {
#		my $val = $val[$i];
#		$val2[$i] = $val;
	#	if ($val =~ /^(8|9)$/ and defined $nopeak{$i}) { # Peak Converted CpG or CH
	#		#$val2[$i] = 7 if $val eq 9;
	#		#$val2[$i] = 6 if $val eq 8;
	#	}
	#}
	#die $print if $totalpeak > 1;
	$print .= "$name\t$totalpeak\n" if $isPeak eq "PEAK";
#	die "$print\n" if $isPeak eq "PEAK";# if $totalpeak == 1;# or $name eq "SEQ_100022";
#	exit 0 if $isPeak eq "PEAK";
	return ($name, \@val, $totalpeak, $peak2{peak});
}

sub find_lots_of_C {
#my ($seqFile, $geneIndex, $box) = @_; #$geneIndexesFa;
#my %geneIndex = %{$geneIndex};
my ($seqFile, $mygene, $outLog) = @_;#, $geneIndex, $box) = @_; #$geneIndexesFa;
my %seq;
print "SEq=$seqFile\n";
print STDERR "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
print $outLog "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
open(my $SEQIN, "<", $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
my $fasta = new FAlite($SEQIN);
my %lotsOfC;

while (my $entry = $fasta->nextEntry()) {
   my $gene = uc($entry->def); $gene =~ s/^>//;
   my $seqz = uc($entry->seq);
	next if $gene ne $mygene;
   print STDERR "\t\tgenez=$gene ($gene)\n";
   print $outLog "\t\tgenez=$gene ($gene)\n";

	my $minlen = 6;
   my $seqz2 = $seqz;#join("", @{$seq{$gene}{seq}});
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
      print "$gene: $beg_C to $end_C ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n";
      $lotsOfC{$gene} .= "C;$beg_C;$end_C,";
   }
	$seqz2 = "";
   $seqz2 =$seqz;# join("", @{$seq{$gene}{seq}});
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
      print "$gene: $beg_G to $end_G ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n";
      $lotsOfC{$gene} .= "G;$beg_G;$end_G,";
   }
	$lotsOfC{$gene} =~ s/,$//;
	print "$lotsOfC{$gene}\n";
	return $lotsOfC{$gene};
}
#foreach my $gene (keys %lotsOfC) {
#   $gene = uc($gene);
#   $lotsOfC{$gene} =~ s/;$//;
#   print "$gene\t$lotsOfC{$gene}\n";
#   my $beg2 = $geneIndex{$gene};
#   foreach my $lines (@{$box->{$gene}}) {
#      print "GENEZ = $gene, lines = $lines\n";
#   }
#   print "genez=$gene,beg=$beg2\n";
#}
#push

}

1;
__END__
# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

   # For nucleotide
# 10 = Nucleotide A
# 11 = Nucleotide C
# 12 = Nucleotide T
# 13 = Nucleotide G






__END__
END1
#my (@inputs, $input1);
#if ($input1 =~ /.orig$/) {
#}
#else {
#	(@inputs) = <$folders/*Pos50.orig>;
#	die "Must have 1 input only! (" . scalar(@inputs) . "):\n" . join("\n-", @inputs) . "\n" if @inputs != 1 and not defined $opt_x;
#	$input1 = defined $opt_x ? $folders : $inputs[0];
#	print "INPUT1=$input1\n";
#}

