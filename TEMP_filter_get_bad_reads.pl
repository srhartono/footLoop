#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_n $opt_S $opt_c $opt_C $opt_o $opt_v $opt_n $opt_b $opt_0);
getopts("s:i:g:f:S:cCo:vn:b:0:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib;
use FAlite;

#my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
#my $footLoopScriptsFolder = dirname(dirname abs_path $0);
#my @version = `$footLoopScriptsFolder/check_software.pl 2>&1 | tail -n 12`;
my ($footLoop_script_folder, $version, $md5script) = check_software();
my $version_small = $version; #"vUNKNOWN";
#my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
	 exit;
}
#my ($version_small) = "vUNKN/CGCGAATGWN";
#foreach my $versionz (@version[0..@version-1]) {
#   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
#}

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N -b $CY<bamFile_bismark_bt2.bam>$N -g $LGN<footLoop_genome.fa>$N -i $LGN<footLoop_geneIndex.bed>$N -o $LGN<output dir>$N

or, if this is a footLoop.pl non-paralel output:
Usage: $YW$0$N -n $CY<folder of -n footLop.pl>$N -o $LGN<output dir>$N

";

(print "\nfootLoop_2_filterBAMFile.pl: $usage\n" and exit 1) unless ex([$opt_b,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1;
#ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1;
(print "\nfootLoop_2_filterBAMFile.pl: please define output (-o)\n" and exit 1) if not defined $opt_o;

my $debug_max_linecount = $opt_0;

my ($footLoop_2_filterBAMFile_outDir) = $opt_o;
makedir($footLoop_2_filterBAMFile_outDir);
my $footLoop_2_filterBAMFile_logFile = "$footLoop_2_filterBAMFile_outDir/footLoop_2_filterBAMFile_logFile.txt";
open (my $outLog, ">", $footLoop_2_filterBAMFile_logFile) or die "Failed to write to footLoop_2_filterBAMFile_logFile: $!\n";

my ($BAMFile, $seqFile, $genez);
if (defined $opt_n) {
	#my ($footLoop_2_filterBAMFile_outDir, $outLog) = @_;
	($BAMFile, $seqFile, $genez) = parse_footLoop_logFile($opt_n, $outLog);
}
else {
	($BAMFile, $seqFile) = ($opt_b, $opt_g);
	$genez = parse_bedFile($opt_i, $outLog);
}

#foreach my $gene (sort keys %{$genez}) {
#	print "$LCY$gene$N: $genez->{$gene}\n";
#}


# check BAM file
LOG($outLog, "\n" . date() . "$YW$0$N Checking files:\n");
LOG($outLog, date() . "Checking BAMFile ($LCY$BAMFile$N)\n");
#check_file($BAMFile, "BAM", $outLog); 

# check seq file
LOG($outLog, date() . "Checking genomeFile ($LCY$seqFile$N)\n");
check_file($seqFile, "seq", $outLog); 

#DIELOG($outLog, "DEBUG Exit before parse_footLoop_logFile\n");
# parse seq file and get chr etc
LOG($outLog, date() . "Parsing genomeFile ($LCY$seqFile$N)\n");
my %refs = %{parse_seqFile($seqFile)}; 
foreach my $chr (sort keys %refs) {
	$genez->{$chr} = @{$refs{$chr}};
}


#DIELOG($outLog, "DEBUG Exit before pasing BAM\n");
my %out;
my %data; my $cons; my %strand;
my $linecount = 0;
my ($BAMFolder, $BAMName) = getFilename($BAMFile, "folderfull");
my $debugFile = "$footLoop_2_filterBAMFile_outDir/debug.txt";
my $outFixedfile = "$footLoop_2_filterBAMFile_outDir/$BAMName.fixed.gz";
my $outBadfile = "$footLoop_2_filterBAMFile_outDir/$BAMName.bad.gz";
my $outWrongMapfile = "$footLoop_2_filterBAMFile_outDir/$BAMName.wrongmap.gz";
open (my $outFixed, "| gzip > $outFixedfile") or die "Cannot write to $outFixedfile: $!\n";
open (my $outBad, "| gzip > $outBadfile") or die "Cannot write to $outBadfile: $!\n";
open (my $outWrongMap, "| gzip > $outWrongMapfile") or die "Cannot write to $outWrongMapfile: $!\n";
open (my $outdebug, ">", "$debugFile") or die "Cannot write to $debugFile: $!\n";
my ($total_read) = `awk '\$2 == 0|| \$2 == 16 {print}' $BAMFile | wc -l` =~ /^\s*(\d+)$/;
$linecount = 0;
my $in1;
if ($BAMFile =~ /.bam$/) {
	open ($in1, "samtools view $BAMFile|") or die "Cannot read from $BAMFile: $!\n" if $BAMFile =~ /.bam$/;
}
else {
	open ($in1, "<", $BAMFile) or die "Cannot read from $BAMFile: $!\n";
}
my $printed = 0;
my %allcount;
LOG($outLog, "\n\n" . date() . "Parsing BAMFile ($LCY$BAMFile$N)\n\n");
my ($maxinsperc, $maxinspercread, $maxdelperc, $maxdelpercread) = (0,"",0,"");
my ($meaninsperc, $meandelperc, $totalread) = (0,0,0);
my %total;
#($total{Read}, $total{Fixed}, $total{Bad}, $total{WrongMap}) = (0,0,0,0);

while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	LOG($outLog, date() . "Done $linecount\n") if $linecount % 100 == 0;
	$linecount ++;
	#last if $linecount == 500;
	#DIELOG($outLog, "DEBUG (-0): Exit at linecount 20\n") 
	last if (defined $opt_0 and $linecount > $debug_max_linecount);

	LOG($outLog, date() . "\t$0: Parsed $LGN$linecount$N / $LCY$total_read$N\n","NA") if $linecount % 50 == 0;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted, @others) = @arr;
	$chr = uc($chr);
	my $others = join("\t", @others); $others = @others == 0 ? "" : "\t$others";
	my @ref1 = defined $refs{$chr} ? @{$refs{$chr}} : die "Can't find gene $chr in $seqFile!\n";
	my @seq1 = split("", $seqs);
	
	# $ref1/$seq1 is original ref/seq from bam line
	# $ref2/$seq2 is $ref1/$seq1 modified by CIGAR
	# $ref3/$seq3 is $ref2/$seq2 modified by bad_region
	my ($ref2, $seq2, $poz, $seqborder0, $seqborder1) = parse_BAMline($line, \@ref1, $outLog);
	
	my %poz = %{$poz};
	my %bad = %{mask_indel_region($ref2, $seq2, $seqborder0, $seqborder1)};
	# mask indel region:
	# GGTAGAG
	# GGT-AGG (G->A)
	# mask this coz the alignment could be GGTAG-G which means there's no G->A anymore
	# instead of being biased we just mask this
	my ($ref3, $seq3, $bad3);
	my $ins = 0;
	my $del = 0;
	my $total;
	my %count;
	my %bigindel;
	my ($indelbeg, $indelend) = (0,0);
	my $isIndel = 0;
	
	my $isBad = 0;
	my $isBadPos = 0;
	my $isWrongMap = 0;
	my $isWrongMapPos = 0;
	my @CT;
	for (my $i = 0; $i < @{$ref2}; $i++) {
		if ($ref2->[$i] ne "-") {
			my $CTtemp;
			if ($seq2->[$i] eq "-") {
				$CTtemp = 0;
			}
			elsif (($ref2->[$i] eq "C" and $seq2->[$i] eq "T") or ($ref2->[$i] eq "G" and $seq2->[$i] eq "A")) {
				$CTtemp = $seq2->[$i];
			}
			elsif ($ref2->[$i] eq "C" or $ref2->[$i] eq "G") {# and $seq2->[$i] eq "T") or ($ref2->[$i] eq "G" and $seq2->[$i] eq "A")) {
					$CTtemp = $ref2->[$i] eq $seq2->[$i] ? $seq2->[$i] : 1;
			}
			else {
				$CTtemp = 1;
			}
			push(@CT, $CTtemp);
		}
		if ($i >= $seqborder0 and $i < $seqborder1) {
			$ins ++ if $ref2->[$i] eq "-";
			$del ++ if $seq2->[$i] eq "-";
			
			# indels
			if ($ref2->[$i] eq "-") {
				if ($isIndel eq 0) {
					$indelbeg = $i;
					$indelend = $i;
					$isIndel = 1;
				}
				else {
					$indelend = $i;
				}
			}
			else {
				if ($isIndel eq 1) {
					$bigindel{$indelbeg} = $indelend;
					$isIndel = 0;
				}
				else {} #do nothing
			}
			
			#if ($ref2->[$i] ne "-" and $seq2->[$i] ne "-") {
			my $ref2nuc = $ref2->[$i];
			my $seq2nuc = $seq2->[$i];
			$count{$ref2nuc}{$seq2nuc} ++;
			#}
			$total ++; # if $ref2->[$i] ne "-";
		}
		my $bad2 = ($i < $seqborder0 or $i >= $seqborder1) ? " " : defined $bad{$i} ? $LRD . "!" . $N : " ";
		if ($ref2->[$i] ne "-") {
			push(@{$ref3}, $ref2->[$i]);
			push(@{$seq3}, $seq2->[$i]);
			push(@{$bad3}, $bad2);
		}
	}
	
	foreach my $tempbeg (sort {$a <=> $b} keys %bigindel) {
		my $tempend = $bigindel{$tempbeg};
		my $templength = $tempend - $tempbeg + 1;
		($isBadPos, $isBad) = ($tempbeg, $templength) if $tempbeg < 400 and $templength > 20;
		($isWrongMapPos, $isWrongMap) = ($tempbeg, $templength) if $tempbeg >= 400 and $templength > 20;
		#print "$read\tins$templength\t$chr:$tempbeg-$tempend\n";
	}
	#my $isBadprint = $isBad ne 0 ? " PromoterBad: pos=$LGN$isBadPos$N len=$LGN${isBad}$N" : " PromoterBad pos=${LGN}0$N len=${LGN}0$N";
	my $isBadprint = " ${LCY}PromoterBad_$LGN$isBadPos$N\_$LGN${isBad}$N";
	my $isWrongMapprint = " ${LPR}WrongMap_$LGN$isWrongMapPos$N\_$LGN${isWrongMap}$N";
#	my $isWrongMapprint = $isWrongMap ne 0 ? " WrongMap_$isWrongMapPos\_${isWrongMap}bp" : "";
	my $insperc = $total == 0 ? 0 : int($ins/$total*10000)/100;
	my $delperc = $total == 0 ? 0 : int($del/$total*10000)/100;
	my $snpperc = "";
	my $indelperc = "";
	my @nuc = qw(- A C G T);
	foreach my $ref2nuc (sort @nuc) {#keys %count) {
		foreach my $seq2nuc (sort @nuc) {#keys %{$count{$ref2nuc}}) {
			my $nuc = $count{$ref2nuc}{$seq2nuc}; $nuc = 0 if not defined $nuc;
			my $perc = $total == 0 ? 0 : int($nuc/$total*10000)/100;
			my $ref2nucprint = $ref2nuc eq "-" ? "del" : $ref2nuc;
			my $seq2nucprint = $seq2nuc eq "-" ? "ins" : $seq2nuc;
			$indelperc .= " $LCY$ref2nucprint$N\_$LPR$seq2nucprint$N=$LGN$perc$N" if $ref2nuc eq "-" or $seq2nuc eq "-";
			$snpperc .= " $LCY$ref2nucprint$N\_$LPR$seq2nucprint$N=$LGN$perc$N" if $ref2nuc ne "-" and $seq2nuc ne "-";
			$allcount{$ref2nuc}{$seq2nuc}{perc} += $perc;
			if (not defined $allcount{$ref2nuc}{$seq2nuc}{maxperc}) {
				$allcount{$ref2nuc}{$seq2nuc}{maxperc} = $perc;
				$allcount{$ref2nuc}{$seq2nuc}{maxread} = $read;
			}
			elsif ($allcount{$ref2nuc}{$seq2nuc}{maxperc} < $perc) {
				$allcount{$ref2nuc}{$seq2nuc}{maxperc} = $perc;
				$allcount{$ref2nuc}{$seq2nuc}{maxread} = $read;
			}
		}
	}
	$meaninsperc += $insperc;
	$meandelperc += $delperc;
	$totalread ++;
	if ($maxinsperc < $insperc) {
		$maxinsperc = $insperc;
		$maxinspercread = $read;
	}
	if ($maxdelperc < $delperc) {
		$maxdelperc = $delperc;
		$maxdelpercread = $read;
	}
	if ($printed < 10) {
		my @ref1print = @ref1 < 10000 ? @ref1 : @ref1[0..10000];
		my @seq1print = @seq1 < 10000 ? @seq1 : @seq1[0..10000];
		my $ref1print = join("", @ref1print);
		my $seq1print = join("", @seq1print);
		#my ($ref1print, $seq1print) = colorconv(\@ref1print, \@seq1print);
		
		LOG($outLog, "Example: $LGN$read$N ($LCY$chr $LGN$pos$N) (ins=$LGN$ins/$total$N ($LGN$insperc$N %), del=$LGN$del/$total$N ($LGN$delperc$N %) $indelperc$snpperc$isBadprint$isWrongMapprint\n\n");
		LOG($outLog, "Original:\n");
		LOG($outLog, "ref: $ref1print\n");
		LOG($outLog, "seq: $seq1print\n\n");

		my @ref2print = @{$ref2} < 10000 ? @{$ref2} : @{$ref2}[0..10000];
		my @seq2print = @{$seq2} < 10000 ? @{$seq2} : @{$seq2}[0..10000];
		my ($ref2print, $seq2print) = colorconv(\@ref2print, \@seq2print);
		LOG($outLog, "After CIGAR:\n");
		LOG($outLog, "ref: $ref2print\n");
		LOG($outLog, "seq: $seq2print\n\n");

		my @ref3print = @{$ref3} < 10000 ? @{$ref3} : @{$ref3}[0..10000];
		my @seq3print = @{$seq3} < 10000 ? @{$seq3} : @{$seq3}[0..10000];
		my @bad3print = @{$bad3} < 10000 ? @{$bad3} : @{$bad3}[0..10000];
		my ($ref3print, $seq3print) = colorconv(\@ref3print, \@seq3print);
		my $bad3print = join("", @bad3print);
		
		LOG($outLog, "After masking those in bad regions:\n");
		LOG($outLog, "ref3: $ref3print\n");
		LOG($outLog, "seq3: $seq3print\n");
		LOG($outLog, "bad3: $bad3print\n\n");
		$printed ++ ;
	}
	
	my	($CTcons, $CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1) = det_C_type($ref3, $seq3, $bad3, $seqborder0, $seqborder1);
	my ($refPrint, $seqPrint) = colorconv($ref3, $seq3);
	my $CTPrint = join("", @{$CTcons});
	my $newstrand = $CT1 > $GA1 ? 0 : $GA1 > $CT1 ? 16 : $strand;

	# 3. DETERMINING TYPE BASED ON CONVERSION
	my $type;

	# 3a 3_NONE: super low C->T conversion and G->A conversion then it's NONE
	if ($CT1 <= 5 and $GA1 <= 5) {
		$type = "3_NONE";
	}

	# 3b. 6_BOTH: otherwise, if C->T and G->A are within +/- 10% then of each other then it's BOTH
	elsif ($CT1 == $GA1 or ($CT1 > 5 and $GA1 > 5 and ($GA1 >= $CT1 * 0.9 and $GA1 <= $CT1 * 1.1) and ($CT1 >= $GA1 * 0.9 and $CT1 <= $GA1 * 1.1))) {
		$type = "6_BOTH";
	}
	# 3c. 1_SNEG and 6_SPOS: otherwise if one is strongly less than the other then STRONG "S" POS or NEG (SPOS or SNEG)
	# -> arbitrary criterias (CT vs GA for POS, and vice versa for NEG)
	#    1. CT 15+ vs GA 5-
	#    2. CT 5+ and ratio CT:GA is at least 3:1
	#    3. CT 20+ and ratio CT:GA is at least 2:1
	elsif (($CT1 > 15 and $GA1 < 5) or ($GA1 >= 5 and $CT1 / $GA1 > 3) or ($GA1 >= 20 and $CT1 / $GA1 >= 2)) {
		$type = "5_SPOS";
	}
	elsif (($GA1 > 15 and $CT1 < 5) or ($CT1 >= 5 and $GA1 / $CT1 > 3) or ($CT1 >= 20 and $GA1 / $CT1 >= 2)) {
		$type = "1_SNEG";
	}
	# 3d. 2_WNEG and 4_WPOS: otherwise just WEAK "W" POS or NEG (WPOS or WNEG)
	elsif ($CT1 < $GA1) {
		$type = "2_WNEG";
	}
	elsif ($CT1 > $GA1) {
		$type = "4_WPOS";
	}
	# 3e. 99_UNK: otherwise default is UNKNOWN
	else {
		$type = "99_UNK";
	}

	$total{$chr}{Read} = 0 if not defined $total{$chr}{Read};
	$total{$chr}{WrongMapBad} = 0 if not defined $total{$chr}{WrongMapBad};
	$total{$chr}{Bad} = 0 if not defined $total{$chr}{Bad};
	$total{$chr}{WrongMap} = 0 if not defined $total{$chr}{WrongMap};
	$total{$chr}{Fixed} = 0 if not defined $total{$chr}{Fixed};

	$total{$chr}{WrongMapBad} ++ if $isBad > 0 and $isWrongMap > 0;

	if ($isBad > 0) {
		$total{$chr}{Bad} ++ if $isWrongMap eq 0;
		print $outBad "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\t" . join("", @CT) . "\n";
	}
	elsif ($isWrongMap > 0) {
		$total{$chr}{WrongMap} ++ if $isBad eq 0;
		print $outWrongMap "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\t" . join("", @CT) . "\n";
	}
	else {
		$total{$chr}{Fixed} ++;
		print $outFixed "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\t" . join("", @CT) . "\n";
	}
	$total{$chr}{Read} ++;
	LOG($outLog, date() . "file=$LCY$BAMFile$N, linecount=$linecount, read=$read, die coz no info\n") if not defined $GG1;

	## 3f. Below is for debug printing
	#print $outdebug ">$read,$type,OldStrand=$strand,NewStrand=$newstrand,$chr,CT0=$CT0,CC0=$CC0,GA0=$GA0,GG0=$GG0,CT1=$CT1,CC1=$CC1,GA1=$GA1,GG1=$GG1\n";
	#print $outdebug "$refPrint\n";
	#print $outdebug "$seqPrint\n";
	#print $outdebug "$CTPrint\n";
	#print $outdebug "$read\t$chr\tstrand=$strand, new=$newstrand\n" . join("", @{$ref3}) . "\n";
	#print $outdebug "CC = $CC1 / $CC0\n";
	#print $outdebug "CT = $CT1 / $CT0\n";
	#print $outdebug "GG = $GG1 / $GG0\n";
	#print $outdebug "GA = $GA1 / $GA0\n";
	#print $outdebug "REF: $refPrint\n";
	#print $outdebug "SEQ: $seqPrint\n";
	#print $outdebug "CON: " . join("", @{$CTcons}) . "\n";
}
close $outBad;
close $outWrongMap;
close $outFixed;
LOG($outLog, "\n");
LOG($outLog, "\n");
foreach my $chr (sort keys %total) {
	LOG($outLog, "chr");
	foreach my $temptype (sort keys %{$total{$chr}}) {
		LOG($outLog, "\t$temptype");
	}
	LOG($outLog, "\n"); last;
}
foreach my $chr (sort keys %total) {
	LOG($outLog, "$chr");
	foreach my $temptype (sort keys %{$total{$chr}}) {
		LOG($outLog, "\t$total{$chr}{$temptype}");
	}
	LOG($outLog, "\n");
}
LOG($outLog, "\n");
LOG($outLog, "\n");


LOG($outLog, "maxinsperc\t$maxinsperc\t$maxinspercread\n");
LOG($outLog, "maxdelperc\t$maxdelperc\t$maxdelpercread\n");
foreach my $ref2nuc (sort keys %allcount) {
	foreach my $seq2nuc (sort keys %{$allcount{$ref2nuc}}) {
		my $maxperc = $allcount{$ref2nuc}{$seq2nuc}{maxperc};
		my $maxread = $allcount{$ref2nuc}{$seq2nuc}{maxread};
		LOG($outLog, "max\t$ref2nuc>$seq2nuc\t$maxperc\t$maxread\n");
	}
}

$meaninsperc = $totalread == 0 ? 0 : int($meaninsperc/$totalread*100)/100;
$meandelperc = $totalread == 0 ? 0 : int($meandelperc/$totalread*100)/100;
LOG($outLog, "meaninsperc\t$meaninsperc\n");
LOG($outLog, "meandelperc\t$meandelperc\n");
foreach my $ref2nuc (sort keys %allcount) {
	foreach my $seq2nuc (sort keys %{$allcount{$ref2nuc}}) {
		my $perc = $allcount{$ref2nuc}{$seq2nuc}{perc};
		$perc = $totalread == 0 ? 0 : int($perc/$totalread * 100)/100;
		LOG($outLog, "mean\t$ref2nuc>$seq2nuc\t$perc\n");
	}
}

#DIELOG($outLog, "DEBUG Exit before printing\n");
foreach my $strand (sort keys %strand) {
	my @types = ("BAMe","diff");
	print $outdebug "$strand: ";
	foreach my $type (@types[0..1]) {
		my $total = $strand{$strand}{$type}; $total = 0 if not defined $total;
		print $outdebug "$type=$total,";
		my $CT = $strand{$strand}{CT}{$type};
		my ($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $CT) {
			$tmm    = int(1000*tmm(@{$CT})+0.5)/1000;
			$mean   = int(1000*mean(@{$CT})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$CT})+0.5)/1000;
			$meanse = int(1000*se(@{$CT})+0.5)/1000;
		}
		print $outdebug "CT=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse, ";
		my $tot = $strand{$strand}{tot}{$type}; 
		($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $tot) {
			$tmm    = int(1000*tmm(@{$tot})+0.5)/1000;
			$mean   = int(1000*mean(@{$tot})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$tot})+0.5)/1000;
			$meanse = int(1000*se(@{$tot})+0.5)/1000;
		}
		print $outdebug "tot=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse\n";
	}
}
exit 0;

# light quick dirty check if BAM or seq file are sane
# BAM is sane if there are at least 10 rows (or less if less than 20 reads) with more than 10 columns
# seq is sane if header is followed by seq and seq is ACTGUN (case-insensitive) at least 20 reads (or less if less than 20 reads in file)
sub check_software {
   my ($footLoop_script_folder, $version, $md5script);
   my @check_software = `check_software.pl 2>&1`;
   foreach my $check_software_line (@check_software[0..@check_software-1]) {
      chomp($check_software_line);
      next if $check_software_line !~ /\=/;
      my ($query, $value) = split("=", $check_software_line);
      next if not defined $query;
      #print "$check_software_line\n";
      if ($query =~ /footLoop_version/) {
         ($version) = $value;
      }
      if ($query =~ /footLoop_script_folder/) {
         next if defined $footLoop_script_folder;
         ($footLoop_script_folder) = $value;
      }
      if ($query =~ /md5sum_script/) {
         ($md5script) = $value;
      }
   }

   print "\ncheck_software.pl\n";
   print "footLoop_script_folder=$footLoop_script_folder\n";
   print "footLoop_version=$version\n";
   print "md5script=$md5script\n\n";
   return($footLoop_script_folder, $version, $md5script);
   #my ($footLoop_script_folder, $version, $md5script) = check_software();
}

sub check_file {
	my ($file, $type, $outLog) = @_;
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: $type file $file does not exist!\n") if ex($file) == 0;
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: $type file $file is empty!\n")       if -s $file  == 0;

	my $filetype = `file -b --mime-type $file`; chomp($filetype);
	my $cmd = ($file =~ /\.(rmdup|bam)$/ or $filetype =~ /(gzip|binary)/) ? "samtools view $file|" : "$file";
	my ($linecount, $check) = (0,0);
	my $currseq; #seq only

	my $total_line = 0;
	my $checkfileIn;
	if ($file =~ /\.(rmdup|bam)$/) {
		($total_line) = `samtools view $file| wc -l` =~ /^\s*(\d+)/;
		LOG($outLog, "samtools view $file| wc -l = $total_line\n","NA");
		open ($checkfileIn, "samtools view $file|") or DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n");
	}
	elsif ($file =~ /\.gz$/ or ($filetype =~ /(gzip|binary)/ and $file !~ /\.(rmdup|bam)$/)) {
		($total_line) = `zcat < $file| wc -l` =~ /^\s*(\d+)/;
		open ($checkfileIn, "zcat < $file|") or DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n");
	}
	else {
		($total_line) = `wc -l $file` =~ /^\s*(\d+)/;
		open ($checkfileIn, "<", $file) or DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n") 
	}
	
	while (my $line = <$checkfileIn>) {
		$linecount ++;
		chomp($line); my @arr = split("\t", $line);
		last if $check >= 20 or $linecount >= $total_line;
		if ($type eq "BAM") {
			if (@arr > 10) {
				$check = $check == 0 ? 2 : $check + 1;
			}
		}
		elsif ($type eq "seq") {
			LOG($outLog, "linecount = $linecount check=$check line = $line\n","NA");
			while ($line !~ /^>/) {
				last if $check >= 20 or $linecount >= $total_line;
				LOG($outLog, "\tlinecount=$linecount, check=$check\n", "NA");
				$currseq .= $line;
				$line = <$checkfileIn>; 
				chomp($line); $linecount ++;
			}
			last if $check >= 20 or $linecount >= $total_line;
			if (defined $currseq and $currseq =~ /^[ACTGUN]+$/i) {
				LOG($outLog, __LINE__ . "before: linecount=$linecount check=$check file=$file type=$type\n","NA");
				$check = $check == 0 ? 0 : $check + 1;
				LOG($outLog, __LINE__ . "after: linecount=$linecount check=$check\n","NA");
				DIELOG($outLog, __LINE__ .  "footLoop_2_filterBAMFile.pl: linecount=$linecount file=$file type=$type, check=$check, check_file failed (does not seem to be a seq file\n\t-> fasta file start with non-header!\n\n") if $check == 0;
				DIELOG($outLog, __LINE__ .  "footLoop_2_filterBAMFile.pl: linecount=$linecount file=$file type=$type check=$check line=$line corruped fasta file!\n\t-> multiple headers in a row\n\n") if $check % 2 != 0;
				undef $currseq;
			}
			if ($linecount < $total_line and defined $line and $line =~ /^>/) {
				LOG($outLog, __LINE__ . "linecount=$linecount check=$check file=$file type=$type\n","NA");
				$check ++;
				# line will always be parsed header (1,3,5,etc) then seq (2,4,6,etc) so if header is even number then corrupted fasta file
				DIELOG($outLog, __LINE__ .  "footLoop_2_filterBAMFile.pl: linecount=$linecount file=$file type=$type check=$check line=$line corruped fasta file!\n") if $check % 2 == 0;
			}
		}
		last if $check >= 20 or $linecount >= $total_line;
	}
	DIELOG($outLog, __LINE__ .  "footLoop_2_filterBAMFile.pl: linecount=$linecount file=$file type=$type, check=$check, check_file failed (does not seem to be a BAM file\n\t-> there's no read with more than 10 collumns!\n") if $check == 0;
}

sub parse_seqFile {
   my ($seqFile) = @_;
   open (my $in, "<", $seqFile) or die "Failed to open seqFile $CY$seqFile$N: $!\n";
	my %ref;
   my $fasta = new FAlite($in);
   while (my $entry = $fasta->nextEntry()) {
      my $def = $entry->def; $def =~ s/>//; $def = uc($def);
      my $seq = uc($entry->seq);
		my @seq = split("", $seq);
		@{$ref{$def}} = @seq;
		LOG($outLog, date() . "REF $def length=" . length($seq) . "\n");# . " SEQ = @seq\n\n");
   }
   close $in;
	return(\%ref);
}

sub parse_footLoop_logFile {
	my ($footLoop_folder, $outLog) = @_;
	$footLoop_folder = "footLoop_folder_unknown" if not defined $opt_n;
#	my $footLoop_folder_forLog = $footLoop_folder;
#	$footLoop_folder_forLog =~ s/\/+/_/g;
#	$footLoop_folder_forLog =~ s/^\/+/SLASH_/;
#my $tempLog = "./.$footLoop_folder_forLog\_TEMPOUTLOG.txt";
#system("touch $tempLog") == 0 or print "Failed to write to $tempLog!\n";
#open (my $tempLogOut, ">", $tempLog) or print "Failed to write to $tempLog: $!\n";
#DIELOG($tempLogOut, "\nTEMPLOG:\n$tempLog\n$footLoop_folder: footLoop_2_filterBAMFile.pl: $usage") unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1 and -e $tempLog;
#DIELOG($tempLogOut, "\nTEMPLOG:\n$tempLog\n$footLoop_folder: footLoop_2_filterBAMFile.pl: please define output (-o)\n") if not defined $opt_o and -e $tempLog;


###########
# LOGFILE #
###########

# log file
my $footLoop_logFile = "$footLoop_folder/logFile.txt";
	LOG($outLog, date() . "Logfile = $LRD$footLoop_logFile$N\n"); 

# parse footLoop logfile

	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: \n\nCan't find $footLoop_logFile! Please run footLoop.pl first before running this!\n\n") if not -e $footLoop_logFile;
#ADD
	$footLoop_logFile = "$footLoop_folder/.PARAMS";
	LOG($outLog, date() . "\t$YW$0$N: LOGFILE=$footLoop_logFile\n");
	my ($BAMFile, $seqFile, $genez);
	if (-e $footLoop_logFile) {
		my @line = `cat $footLoop_logFile`;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i]; chomp($line);
			if ($line =~ /footLoop.pl,BAMFile/) {
				($BAMFile) = $line =~ /BAMFile,(.+)$/;
			}
			if ($line =~ /footLoop.pl,seqFile/) {
				($seqFile) = $line =~ /seqFile,(.+)$/;
			}
			if ($line =~ /footLoop.pl,geneIndexFile/) {
				my ($geneIndexFile) = $line =~ /geneIndexFile,(.+)$/;
				my @line = `cut -f1-4 $geneIndexFile`;
				foreach my $line (@line) {
					chomp($line); next if $line =~ /^#/;
					my ($chr, $beg, $end, $gene) = split("\t", $line);
					my $length = $end - $beg;
					$gene = uc($gene);
					$genez->{$gene} = $length;
					LOG($outLog, "\t$YW\$0${N}::${LGN}parse_footLoop_logFile$N: gene $gene = $length bp\n");
				}
			}
			
#			if ($line =~ /^!BAMFile=/) {
#				($BAMFile) = $line =~ /^!BAMFile=(.+)$/;
#			}
#			if ($line =~ /^!seqFile=/) {
#				($seqFile) = $line =~ /^!seqFile=(.+)$/;
#			}
#			if ($line =~ /gene=.+length=/) {
#				my ($gene, $length) = $line =~ /^.+gene=(.+) length=(\d+)$/;
#				die if not defined $gene or not defined $length;
#				$genez->{$gene} = $length;
#				print "gene $gene = $length bp\n";
#			}
#			last if $line =~ /footLoop_2_filterBAMFile.pl/;
		}
	}
	return($BAMFile, $seqFile, $genez);
}


sub parse_bedFile {
	my ($input1, $outLog) = @_;
	LOG($outLog, date() . "$YW${0}${N}::parse_bedFile: Parsing bedFile from $LCY$input1$N\n");
	my $genez;
	open (my $in1, "<", $input1) or DIELOG($outLog, "Failed to read from $input1: $!\n");
	while (my $line = <$in1>) {
		chomp($line); next if $line =~ /^#/;
		my ($chr, $beg, $end, $gene) = split("\t", $line);
		my $length = $end - $beg;
		$gene = uc($gene);
		$genez->{$gene} = $length;
		LOG($outLog, "\t$YW\$0${N}::${LGN}parse_bedFile$N: gene $gene = $length bp\n");
	}
	close $in1;
	return($genez);
}


sub mask_indel_region {

	#CCAGA
	#CC-GA
	#make C-GA position as "bad"
	# because the deletion might mess with the alignment
	my ($ref2, $seq2, $seqborder0, $seqborder1) = @_;
	my %bad;
	my @bad; my $prints;
	for (my $i = $seqborder0; $i < $seqborder1; $i++) {
		my $badcount = 0;
		my $reftemp = "";
		my $seqtemp = "";
		for (my $j = $i; $j < $seqborder1; $j++) {
			DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Undefined ref2 at j=$j\n") if not defined $ref2->[$j];
			$reftemp .= $ref2->[$j];
			DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Undefined seq2 at j=$j\n") if not defined $seq2->[$j];
			$seqtemp .= $seq2->[$j];
			if ($ref2->[$j] eq "-" or $seq2->[$j] eq "-") {
				$badcount ++;
			}
			else {
				last;
			}
		}
		my $beg = $i - $badcount < 0 ? 0 : $i - $badcount;
		my $end = $i + 2*$badcount > @{$ref2} ? @{$ref2} : $i + 2*$badcount;
		$prints .= join("", (" ")x$beg) . join("", ("!")x($end-$beg+1)) . "($beg-$end)\n" if $badcount > 0;
		$prints .= join("", (" ")x$i) . "$reftemp\n" if $badcount > 0;
		$prints .= join("", (" ")x$i) . "$seqtemp\n" if $badcount > 0;
		for (my $j = $beg; $j <= $end; $j++) {
			last if $badcount == 0;
			$bad[$j] = 0;
			$bad{$j} = 1;
		}
		$i += $badcount;
	}
	for (my $i = 0; $i < @bad; $i++) {
		$bad[$i] = " " if not defined $bad[$i];
	}
	return(\%bad);
}

sub parse_BAMline {
	my ($line, $refs, $outLog) = @_;
	my @refs = @{$refs};
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted) = split("\t", $line);
	my @seq = split("",$seqs);
	my ($num, $alp, $lengthseq) = parse_cigar($cigar); die if not defined $num;
	my @num  = @{$num}; my @alp = @{$alp};
	my @ref0 = @refs[0..$pos-2];
	my @ref = @refs[$pos-1..@refs-1];
	my ($seq, $ref, $seqpos, $refpos) = (\@seq, \@ref, 0, 0);
	my %pos;
	my $lengthref = @ref; my $insref = 0; my $insseq = 0;
	for (my $i = 0; $i < @num; $i++) {
		($ref) = $alp[$i] eq "I" ? ins($ref, $refpos, "-", $num[$i], "ref") : $ref;
		($seq) = $alp[$i] eq "D" ? ins($seq, $seqpos, "-", $num[$i], "seq") : $seq; 
		$refpos += $num[$i];
		$seqpos += $num[$i];
		$insref += $num[$i] if $alp[$i] eq "I";
		$insseq += $num[$i] if $alp[$i] eq "D";
	}
	my $refend = $refpos - $insref;
	@ref = (@ref0, @{$ref});
	my $seqborder0 = $pos-1;
	my $seqborder1 = $pos - 1 + @{$seq};
	@seq = ((("-") x ($pos-1)), @{$seq}, (("-")x($lengthref-$refend)));
	return(\@ref, \@seq, \%pos, $seqborder0, $seqborder1);
}

sub det_C_type {
   my ($ref, $seq, $bad, $seqborder0, $seqborder1) = @_;
   my @ref = @{$ref};
   my @seq = @{$seq};
	my @bad = @{$bad};
   my (@top);
   @ref = ("-","-",@ref,"-","-");
   @seq = ("-","-",@seq,"-","-");
   @bad = (" "," ",@bad," "," ");
   $seqborder0 += 2;
   $seqborder1 += 2;
   my $len = @ref;
	my ($CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1) = (0,0,0,0,0,0,0,0);
   for (my $i = 2; $i < $len-2; $i++) {
      my ($beg, $end, $pos) = ($i-2, $i+2, $i-2); #pos is for top and bot position
      my ($top1, $top2, $top3, $top4, $top5) = @ref[$beg..$end]; #NNCNN
      my ($bot1, $bot2, $bot3, $bot4, $bot5) = @seq[$beg..$end]; #NNGNN
		my $chunkbot = join("", @seq[$beg..$end]);
      if ($i < $seqborder0 or $i >= $seqborder1) {
			$top[$pos] = "-";
      }
      elsif ($top3 eq "C") {
			my $top;
			if ($bad[$i] eq " ") {
	         $top = $top4 eq "G" ? "E" : $top5 eq "G" ? "D" : "C";
         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
			else {
	         $top = $top4 eq "G" ? "3" : $top5 eq "G" ? "2" : "1";
         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
         $top = $bot3 eq "T" ? lc($top) : ($top =~ /^(C|D)$/ and $bot3 eq "-") ? "B" : ($top =~ /^E$/ and $bot3 eq "-") ? "F" : $top;
         $top =~ tr/CDE123/MNOPQR/ if $bot3 !~ /^(C|T)$/; #not CC or CT.
			#CH=B(-)CDcdMN, CG=F(-)EeO
         $top[$pos] = $top;
      }
      elsif ($top3 eq "G") {
			my $top;
			if ($bad[$i] eq " ") {
	         $top = $top2 eq "C" ? "I" : $top1 eq "C" ? "H" : "G";
	         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
			else {
	         $top = $top2 eq "C" ? "6" : $top1 eq "C" ? "5" : "4";
	         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
         $top = $bot3 eq "A" ? lc($top) : ($top =~ /^(G|H)$/ and $bot3 eq "-") ? "J" : ($top eq "I" and $bot3 eq "-") ? "K" : $top;
         $top =~ tr/GHI456/UVWXYZ/ if $bot3 !~ /^(G|A)$/; # not GG or GA. 
			#GH=J(-)GHghUV, GC=K(-)IiW
			$top[$pos] = $top;
      }
      else {
			$top[$pos] = $bot3 eq "-" ? "_" : "$bot3";
      }
		$CC0 ++ if "$top3$bot3" eq "CC"; 
		$GG0 ++ if "$top3$bot3" eq "GG"; 
		$CT0 ++ if "$top3$bot3" eq "CT"; 
		$GA0 ++ if "$top3$bot3" eq "GA";
		if ($bad[$i] eq " ") {
			my $left = $i - 5 < 0 ? 0 : $i - 5;
			my $rite = $i + 5 >= @ref ? @ref - 1: $i + 5;
			my $chunk = join("", @seq[$left..$rite]);
			if ($chunk !~ /\-/) {
				$CT1 ++ if "$top3$bot3" eq "CT";
				$GA1 ++ if "$top3$bot3" eq "GA";
				$CC1 ++ if "$top3$bot3" eq "CC"; 
				$GG1 ++ if "$top3$bot3" eq "GG"; 
			}
		}
   }
   return(\@top, $CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1);
}
sub ins {
	my ($arr, $pos, $ins, $total, $type) = @_;
	my @ins = ($ins) x $total;
	my @arr = @{$arr};
	my @arr0 = @arr[0..$pos-1];
	my @arr1 = @arr[$pos..@arr-1];
	die if @arr0 == 0;
	@arr = (@arr0, @ins, @arr1);
	return(\@arr);
}

sub getConv {
	my ($converted) = @_;
	($converted) = $converted =~ /^.+:([\.A-Z]+)$/i; my $length2 = length($converted);
	my ($X) = $converted =~ tr/X/X/;
	my ($U) = $converted =~ tr/U/U/;
	my ($H) = $converted =~ tr/H/H/;
	my ($Z) = $converted =~ tr/Z/Z/;
	my ($x) = $converted =~ tr/x/x/;
	my ($u) = $converted =~ tr/u/u/;
	my ($h) = $converted =~ tr/h/h/;
	my ($z) = $converted =~ tr/z/z/;
	my ($dot) = $converted =~ tr/\./\./;
	my $length = ($X+$H+$Z+$U+$x+$h+$z+$u+$dot);
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: X=$X, H=$H, Z=$Z, U=$u, x=$x, h=$h, z=$z, u=$u, dot=$dot, length=$length, length2=$length2\n\n") if $length ne $length2;
	$data{stat} = "X=$X, H=$H, Z=$Z, U=$u, x=$x, h=$h, z=$z, u=$u, dot=$dot, lengthsum=$length, lengthcol14=$length2";
	my $conv = $x + $h + $z + $u;
	my $notc = $X + $H + $Z + $U;
	my $nonc = $dot;
	my @converted = split("", $converted);
	return(\@converted);

}
sub check_chr_in_BAM {
	my ($BAMFile) = @_;
	open (my $in2, "cut -f2,3 $BAMFile|") or die "Cannot read from $BAMFile: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		my @arr = split("\t", $line);
		next if @arr == 0;
		next if $arr[0] !~ /^\d+$/;
		$linecount ++;
		my ($strand, $chr) = @arr; $chr = uc($chr);
		next if $strand eq 4;
		DIELOG($outLog, "footLoop_2_filterBAMFile.pl: Can't find gene $chr in $seqFile!\n") if not defined $refs{$chr};
		next;
	}
}
__END__

