#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_n $opt_S $opt_c $opt_C $opt_o $opt_v $opt_n $opt_b $opt_0);
getopts("s:i:g:f:S:cCo:vn:b:0");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib;
use FAlite;

my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);
my $version = "v3.8";
my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $version_small = $version; 
if (defined $opt_v) {
   print "$version\n";
	 exit;
}

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N -b $CY<bamFile_bismark_bt2.bam>$N -g $LGN<footLoop_genome.fa>$N -i $LGN<footLoop_geneIndex.bed>$N -o $LGN<output dir>$N

or, if this is a footLoop.pl non-paralel output:
Usage: $YW$0$N -n $CY<folder of -n footLop.pl>$N -o $LGN<output dir>$N

";

(print "\nfootLoop_2_filterBAMFile.pl: $usage\n" and exit 1) unless ex([$opt_b,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1;
(print "\nfootLoop_2_filterBAMFile.pl: please define output (-o)\n" and exit 1) if not defined $opt_o;

my ($footLoop_2_filterBAMFile_outDir) = $opt_o;
makedir($footLoop_2_filterBAMFile_outDir);
my $footLoop_2_filterBAMFile_logFile = "$footLoop_2_filterBAMFile_outDir/footLoop_2_filterBAMFile_logFile.txt";
$footLoop_2_filterBAMFile_logFile = "$footLoop_2_filterBAMFile_outDir/footLoop_2_filterBAMFile_logFile_debug.txt" if defined $opt_0;
open (my $outLog, ">", $footLoop_2_filterBAMFile_logFile) or die "Failed to write to footLoop_2_filterBAMFile_logFile: $!\n";

my ($BAMFile, $seqFile, $genez);
if (defined $opt_n) {
	($BAMFile, $seqFile, $genez) = parse_footLoop_logFile($opt_n, $outLog);
}
else {
	($BAMFile, $seqFile) = ($opt_b, $opt_g);
	$genez = parse_bedFile($opt_i, $outLog);
}

# check BAM file
LOG($outLog, "\n" . date() . "$YW$0$N Checking files:\n");
LOG($outLog, date() . "Checking BAMFile ($LCY$BAMFile$N)\n");
check_file($BAMFile, "BAM", $outLog); 

# check seq file
LOG($outLog, date() . "Checking genomeFile ($LCY$seqFile$N)\n");
check_file($seqFile, "seq", $outLog); 

# parse seq file and get chr etc
LOG($outLog, date() . "Parsing genomeFile ($LCY$seqFile$N)\n");
my %refs = %{parse_seqFile($seqFile)}; 
foreach my $chr (sort keys %refs) {
	$genez->{$chr} = @{$refs{$chr}};
}


my %out;
my %data; my $cons; my %strand;
my $linecount = 0;
my ($BAMFolder, $BAMName) = getFilename($BAMFile, "folderfull");
my $outFixedfile = "$footLoop_2_filterBAMFile_outDir/$BAMName.fixed.gz";
my $outFixed;
my $outdebug;
if (not defined $opt_0) {
	open ($outFixed, "| gzip > $outFixedfile") or die "Cannot write to $outFixedfile: $!\n";
}
my ($total_read) = 0;
$linecount = 0;
my $in1;
if ($BAMFile =~ /.bam$/) {
	open ($in1, "samtools view $BAMFile|") or die "Cannot read from $BAMFile: $!\n" if $BAMFile =~ /.bam$/;
	($total_read) = `samtools view $BAMFile|wc -l`;
}
else {
	open ($in1, "<", $BAMFile) or die "Cannot read from $BAMFile: $!\n";
	($total_read) = `cat $BAMFile|wc -l`;
}
my $printed = 0;
my %allcount;
LOG($outLog, "\n\n" . date() . "Parsing BAMFile ($LCY$BAMFile$N)\n\n");
my ($maxinsperc, $maxinspercread, $maxdelperc, $maxdelpercread) = (0,"",0,"");
my ($meaninsperc, $meandelperc, $totalread) = (0,0,0);
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	$linecount ++;
	DIELOG($outLog, "DEBUG (-0): Exit at linecount 20\n") if defined $opt_0 and $linecount > 1000;

	LOG($outLog, date() . "\t$0: Parsed $LGN$linecount$N / $LCY$total_read$N\n") if $linecount % 100 == 0;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted, @others) = @arr;
	
	LOG($outLog, "$read,cigar1,$cigar\n","NA");

	$chr = uc($chr);
	my $others = join("\t", @others); $others = @others == 0 ? "" : "\t$others";
	my @ref1 = defined $refs{$chr} ? @{$refs{$chr}} : die "Can't find gene $chr in $seqFile!\n";
	my @seq1 = split("", $seqs);
	
	# $ref1/$seq1 is original ref/seq from bam line
	# $ref2/$seq2 is $ref1/$seq1 modified by CIGAR
	# $ref3/$seq3 is $ref2/$seq2 modified by bad_region
	my ($ref2, $seq2, $poz, $seqborder0, $seqborder1) = parse_BAMline($line, \@ref1, $outLog);
	
	my %poz = %{$poz};
	my %bad;

	# mask indel region:
	# GGTAGAG
	# GGT-AGG (G->A)
	# mask this coz the alignment could be GGTAG-G which means there's no G->A anymore
	# instead of being biased we just mask this
	my ($bad2, $ref3, $seq3, $bad3);
	my $ins = 0;
	my $del = 0;
	my $total;
	my %count;
	for (my $i = 0; $i < @{$ref2}; $i++) {
		if ($i >= $seqborder0 and $i < $seqborder1) {
			$ins ++ if $ref2->[$i] eq "-";
			$del ++ if $seq2->[$i] eq "-";
			my $ref2nuc = $ref2->[$i];
			my $seq2nuc = $seq2->[$i];
			$count{$ref2nuc}{$seq2nuc} ++;
			$total ++; 
		}
		my $badnuc2 = ($i < $seqborder0 or $i >= $seqborder1) ? " " : defined $bad{$i} ? "!" : " ";
		push(@{$bad2}, $badnuc2);
		if ($ref2->[$i] ne "-") {
			push(@{$ref3}, $ref2->[$i]);
			push(@{$seq3}, $seq2->[$i]);
			push(@{$bad3}, $badnuc2);
		}
	}
	
	my $insperc = $total == 0 ? 0 : int($ins/$total*10000)/100;
	my $delperc = $total == 0 ? 0 : int($del/$total*10000)/100;
	foreach my $ref2nuc (sort keys %count) {
		foreach my $seq2nuc (sort keys %{$count{$ref2nuc}}) {
			my $nuc = $count{$ref2nuc}{$seq2nuc};
			my $perc = $total == 0 ? 0 : int($nuc/$total*10000)/100;
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
	if ($printed <= 1000) {
		my @ref1print = @ref1 < 10000 ? @ref1 : @ref1[0..10000];
		my @seq1print = @seq1 < 10000 ? @seq1 : @seq1[0..10000];
		my ($ref1print,$seq1print) = (join("", @ref1print),join("", @seq1print));
		
		my $NA = $printed < 5 ? "" : "NA";
		my $toprint = "";
		$toprint .= "$read,GOOD,$linecount,ref,seq1,$chr,$ref1print\n";
		$toprint .= "$read,GOOD,$linecount,que,seq1,$chr,$seq1print\n";

		my @ref2print = @{$ref2} < 10000 ? @{$ref2} : @{$ref2}[0..10000];
		my @seq2print = @{$seq2} < 10000 ? @{$seq2} : @{$seq2}[0..10000];
		my @bad2print = @{$bad2} < 10000 ? @{$bad2} : @{$bad2}[0..10000];
		my ($ref2print,$seq2print) = (join("", @ref2print),join("", @seq2print));
		my $bad2print = join("", @bad2print);
		$toprint .= "$read,GOOD,$linecount,ref,seq2,$chr,$ref2print\n";
		$toprint .= "$read,GOOD,$linecount,que,seq2,$chr,$seq2print\n";
		$toprint .= "$read,GOOD,$linecount,bad,seq2,$chr,$bad2print\n";

		my @ref3print = @{$ref3} < 10000 ? @{$ref3} : @{$ref3}[0..10000];
		my @seq3print = @{$seq3} < 10000 ? @{$seq3} : @{$seq3}[0..10000];
		my @bad3print = @{$bad3} < 10000 ? @{$bad3} : @{$bad3}[0..10000];
		my ($ref3print,$seq3print) = (join("", @ref3print),join("", @seq3print));
		my $bad3print = join("", @bad3print);
		
		$toprint .= "$read,GOOD,$linecount,ref,seq3,$chr,$ref3print\n";
		$toprint .= "$read,GOOD,$linecount,que,seq3,$chr,$seq3print\n";
		$toprint .= "$read,GOOD,$linecount,bad,seq3,$chr,$bad3print\n";
		$printed ++ ;


		my ($VRmetpos, $VRmetneg, $VRmatperc, $VRmisperc, $VRinsperc, $VRdelperc, $VRseqlength);
		($VRmetpos, $VRmetneg, $VRmatperc, $VRmisperc, $VRinsperc, $VRdelperc, $VRseqlength, $toprint)	 = check_VR($read, $linecount, \@{$ref2}, \@{$seq2}, \@{$ref3}, \@{$seq3}, $outLog, $NA, $toprint);
		my $VRinstotal = int($VRinsperc/100 * $VRseqlength+0.5);
		my $VRdeltotal = int($VRdelperc/100 * $VRseqlength+0.5);
		$toprint = "$read,GOOD,$linecount,all,info,$chr,$pos,$ins,$del,$total,$insperc,$delperc\n" . $toprint . "\n";
		if ($toprint =~ /,big,/) {
			$toprint =~ s/,GOOD,/,BIGINDEL,/g;
		}
		DIELOG($outLog, "GOOD DEBUG\n") if $printed > 10 and defined $opt_0;
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

	if (not defined $opt_0) {
		print $outFixed "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\n";
	}
	LOG($outLog, date() . "file=$LCY$BAMFile$N, linecount=$linecount, read=$read, die coz no info\n") if not defined $GG1;

}
if (not defined $opt_0) {
	close $outFixed;
}

LOG($outLog, "maxinsperc\t$maxinsperc\t$maxinspercread\n","NA");
LOG($outLog, "maxdelperc\t$maxdelperc\t$maxdelpercread\n","NA");
foreach my $ref2nuc (sort keys %allcount) {
	foreach my $seq2nuc (sort keys %{$allcount{$ref2nuc}}) {
		my $maxperc = $allcount{$ref2nuc}{$seq2nuc}{maxperc};
		my $maxread = $allcount{$ref2nuc}{$seq2nuc}{maxread};
		LOG($outLog, "max\t$ref2nuc>$seq2nuc\t$maxperc\t$maxread\n","NA");
	}
}

$meaninsperc = $totalread == 0 ? 0 : int($meaninsperc/$totalread*100)/100;
$meandelperc = $totalread == 0 ? 0 : int($meandelperc/$totalread*100)/100;
LOG($outLog, "meaninsperc\t$meaninsperc\n","NA");
LOG($outLog, "meandelperc\t$meandelperc\n","NA");
foreach my $ref2nuc (sort keys %allcount) {
	foreach my $seq2nuc (sort keys %{$allcount{$ref2nuc}}) {
		my $perc = $allcount{$ref2nuc}{$seq2nuc}{perc};
		$perc = $totalread == 0 ? 0 : int($perc/$totalread * 100)/100;
		LOG($outLog, "mean\t$ref2nuc>$seq2nuc\t$perc\n","NA");
	}
}

foreach my $strand (sort keys %strand) {
	my @types = ("BAMe","diff");
	if (not defined $opt_0) {
		print $outdebug "$strand: ";
	}
	foreach my $type (@types[0..1]) {
		my $total = $strand{$strand}{$type}; $total = 0 if not defined $total;
		if (not defined $opt_0) {
			print $outdebug "$type=$total,";
		}
		my $CT = $strand{$strand}{CT}{$type};
		my ($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $CT) {
			$tmm    = int(1000*tmm(@{$CT})+0.5)/1000;
			$mean   = int(1000*mean(@{$CT})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$CT})+0.5)/1000;
			$meanse = int(1000*se(@{$CT})+0.5)/1000;
		}
		if (not defined $opt_0) {
			print $outdebug "CT=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse, ";
		}
		my $tot = $strand{$strand}{tot}{$type}; 
		($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $tot) {
			$tmm    = int(1000*tmm(@{$tot})+0.5)/1000;
			$mean   = int(1000*mean(@{$tot})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$tot})+0.5)/1000;
			$meanse = int(1000*se(@{$tot})+0.5)/1000;
		}
		if (not defined $opt_0) {
			print $outdebug "tot=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse\n";
		}
	}
}
exit 0;

# light quick dirty check if BAM or seq file are sane
# BAM is sane if there are at least 10 rows (or less if less than 20 reads) with more than 10 columns
# seq is sane if header is followed by seq and seq is ACTGUN (case-insensitive) at least 20 reads (or less if less than 20 reads in file)

sub check_VR {
	my ($read, $linecount, $refs2, $ques2, $refs3, $ques3, $outLog, $NA, $toprint) = @_;
	
	my @return;
	my @ref2 = @{$refs2};
	my @que2 = @{$ques2};
	my $ref2 = join("", @ref2);
	my $que2 = join("", @que2);
	my %dict;
	my @ref3 = @{$refs3};
	my @que3 = @{$ques3};
	my $ref3 = join("", @ref3);
	my $que3 = join("", @que3);
	my $ref2length = length($ref2);
	my $ref3length = length($ref3);
	my ($VRrefseq0, $VRrefseq1) = ("CTCCTCGCCTCGGTCACTGCGACGAATTC","CCGAGGGAGTGTTTGGGTACCCATTATAC");
	my ($VRrefseq0rev, $VRrefseq1rev) = (revcomp($VRrefseq1),revcomp($VRrefseq0));
	my ($leftrefseq, $VRrefseq3, $rightrefseq) = ("","","");
	my ($leftqueseq, $VRqueseq3, $rightqueseq) = ("","","");
	my ($leftrefseq2, $VRrefseq2, $rightrefseq2) = ("","","");
	my ($leftqueseq2, $VRqueseq2, $rightqueseq2) = ("","","");
	my $ref3Ind = 0;
	# ref2: ACTG---AAA
	#       0123456789
	# ref3: ACTG---AAA
	#       0123333456
	#2i  3i
	#3 = 3
	#4 = 3
	#5 = 3
	#6 = 3
	#4 = 7
	#5 = 8
	#
	for (my $ref2Ind = 0; $ref2Ind < @ref2; $ref2Ind++) {
		if (not defined $dict{ref3toref2}{$ref3Ind}) {
			$dict{ref3toref2}{$ref3Ind} = $ref2Ind;
		}
		$dict{ref2toref3}{$ref2Ind} = $ref3Ind;
		$ref3Ind ++ if $ref2[$ref2Ind] ne "-";
	}
	$dict{ref3toref2}{length($ref3)} = length($ref2);
	$dict{ref2toref3}{length($ref2)} = length($ref3);

	# PRIMER
	if (length($ref2) > 200) {
		my $FWrefseq3pos1 = 0;
		my $FWrefseq3pos2 = 100;
		my $RVrefseq3pos1 = length($ref3)-100;
		my $RVrefseq3pos2 = length($ref3);
		my $FWrefseq2pos1 = $dict{ref3toref2}{$FWrefseq3pos1};
		my $FWrefseq2pos2 = $dict{ref3toref2}{$FWrefseq3pos2};
		my $RVrefseq2pos1 = $dict{ref3toref2}{$RVrefseq3pos1};
		my $RVrefseq2pos2 = defined $dict{ref3toref2}{$RVrefseq3pos2} ? $dict{ref3toref2}{$RVrefseq3pos2} : length($ref2);

		my ($FWprimerrefseq2, $RVprimerrefseq2) = (substr($ref2,$FWrefseq2pos1,100),substr($ref2,$RVrefseq2pos1,$RVrefseq2pos2-$RVrefseq2pos1));
		my ($FWprimerqueseq2, $RVprimerqueseq2) = (substr($que2,$FWrefseq2pos1,100),substr($que2,$RVrefseq2pos1,$RVrefseq2pos2-$RVrefseq2pos1));
		my ($FWprimerbadseq2, $FWprimerins2, $FWprimerdel2, $FWprimermat2, $FWprimermis2, $FWprimermetpos2, $FWprimermetneg2, $FWprimerseqlength2) = get_seq_indel_stat($FWprimerrefseq2, $FWprimerqueseq2);
		my $FWprimermatperc2 = int(1000*$FWprimermat2 / length($FWprimerrefseq2)+0.5)/10;
		my $FWprimermisperc2 = int(1000*$FWprimermis2 / length($FWprimerrefseq2)+0.5)/10;
		my $FWprimerinsperc2 = int(1000*$FWprimerins2 / length($FWprimerrefseq2)+0.5)/10;
		my $FWprimerdelperc2 = int(1000*$FWprimerdel2 / length($FWprimerrefseq2)+0.5)/10;
		my ($RVprimerbadseq2, $RVprimerins2, $RVprimerdel2, $RVprimermat2, $RVprimermis2, $RVprimermetpos2, $RVprimermetneg2, $RVprimerseqlength2) = get_seq_indel_stat($RVprimerrefseq2, $RVprimerqueseq2);
		my $RVprimermatperc2 = int(1000*$RVprimermat2 / length($RVprimerrefseq2)+0.5)/10;
		my $RVprimermisperc2 = int(1000*$RVprimermis2 / length($RVprimerrefseq2)+0.5)/10;
		my $RVprimerinsperc2 = int(1000*$RVprimerins2 / length($RVprimerrefseq2)+0.5)/10;
		my $RVprimerdelperc2 = int(1000*$RVprimerdel2 / length($RVprimerrefseq2)+0.5)/10;

		my ($FWprimerrefseq3, $RVprimerrefseq3) = (substr($ref3,$FWrefseq3pos1,100),substr($ref3,$RVrefseq3pos1,$RVrefseq3pos2-$RVrefseq3pos1));
		my ($FWprimerqueseq3, $RVprimerqueseq3) = (substr($que3,$FWrefseq3pos1,100),substr($que3,$RVrefseq3pos1,$RVrefseq3pos2-$RVrefseq3pos1));
		my ($FWprimerbadseq3, $FWprimerins3, $FWprimerdel3, $FWprimermat3, $FWprimermis3, $FWprimermetpos3, $FWprimermetneg3, $FWprimerseqlength3) = get_seq_indel_stat($FWprimerrefseq3, $FWprimerqueseq3);
		my $FWprimermatperc3 = int(1000*$FWprimermat3 / length($FWprimerrefseq3)+0.5)/10;
		my $FWprimermisperc3 = int(1000*$FWprimermis3 / length($FWprimerrefseq3)+0.5)/10;
		my $FWprimerinsperc3 = int(1000*$FWprimerins3 / length($FWprimerrefseq3)+0.5)/10;
		my $FWprimerdelperc3 = int(1000*$FWprimerdel3 / length($FWprimerrefseq3)+0.5)/10;
		my ($RVprimerbadseq3, $RVprimerins3, $RVprimerdel3, $RVprimermat3, $RVprimermis3, $RVprimermetpos3, $RVprimermetneg3, $RVprimerseqlength3) = get_seq_indel_stat($RVprimerrefseq3, $RVprimerqueseq3);
		my $RVprimermatperc3 = int(1000*$RVprimermat3 / length($RVprimerrefseq3)+0.5)/10;
		my $RVprimermisperc3 = int(1000*$RVprimermis3 / length($RVprimerrefseq3)+0.5)/10;
		my $RVprimerinsperc3 = int(1000*$RVprimerins3 / length($RVprimerrefseq3)+0.5)/10;
		my $RVprimerdelperc3 = int(1000*$RVprimerdel3 / length($RVprimerrefseq3)+0.5)/10;
		
		my $FWprimerrefseqprint2 = $FWprimerrefseq2;
		my $RVprimerrefseqprint2 = $RVprimerrefseq2;
		my $FWprimerrefseqprint3 = $FWprimerrefseq3;
		my $RVprimerrefseqprint3 = $RVprimerrefseq3;

		my $FWprimerqueseqprint2 = $FWprimerqueseq2;
		my $RVprimerqueseqprint2 = $RVprimerqueseq2;
		my $FWprimerqueseqprint3 = $FWprimerqueseq3;
		my $RVprimerqueseqprint3 = $RVprimerqueseq3;
	
		$toprint .= "$read,GOOD,$linecount,all,FWprimerpos2,$FWrefseq2pos1,$FWrefseq2pos2,$FWprimermetpos2,$FWprimermetneg2,$FWprimermatperc2,$FWprimermisperc2,$FWprimerinsperc2,$FWprimerdelperc2,$FWprimerseqlength2\n";
		$toprint .= "$read,GOOD,$linecount,ref,FWprimerseq2,$FWprimerrefseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,que,FWprimerseq2,$FWprimerqueseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,bad,FWprimerseq2,$FWprimerbadseq2\n";

		$toprint .= "$read,GOOD,$linecount,all,RVprimerpos2,$RVrefseq2pos1,$RVrefseq2pos2,$RVprimermetpos2,$RVprimermetneg2,$RVprimermatperc2,$RVprimermisperc2,$RVprimerinsperc2,$RVprimerdelperc2,$RVprimerseqlength2\n";
		$toprint .= "$read,GOOD,$linecount,ref,RVprimerseq2,$RVprimerrefseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,que,RVprimerseq2,$RVprimerqueseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,bad,RVprimerseq2,$RVprimerbadseq2\n";

		$toprint .= "$read,GOOD,$linecount,all,FWprimerpos3,$FWrefseq3pos1,$FWrefseq3pos2,$FWprimermetpos3,$FWprimermetneg3,$FWprimermatperc3,$FWprimermisperc3,$FWprimerinsperc3,$FWprimerdelperc3,$FWprimerseqlength3\n";
		$toprint .= "$read,GOOD,$linecount,ref,FWprimerseq3,$FWprimerrefseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,que,FWprimerseq3,$FWprimerqueseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,bad,FWprimerseq3,$FWprimerbadseq3\n";

		$toprint .= "$read,GOOD,$linecount,all,RVprimerpos3,$RVrefseq3pos1,$RVrefseq3pos2,$RVprimermetpos3,$RVprimermetneg3,$RVprimermatperc3,$RVprimermisperc3,$RVprimerinsperc3,$RVprimerdelperc3,$RVprimerseqlength3\n";
		$toprint .= "$read,GOOD,$linecount,ref,RVprimerseq3,$RVprimerrefseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,que,RVprimerseq3,$RVprimerqueseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,bad,RVprimerseq3,$RVprimerbadseq3\n";


	}

	if ($ref3 =~ /$VRrefseq0rev.+$VRrefseq1rev/) {
		($VRrefseq0, $VRrefseq1) = ($VRrefseq0rev, $VRrefseq1rev);
	}

	if ($ref3 =~ /$VRrefseq0.+$VRrefseq1/) {
		($leftrefseq, $VRrefseq3, $rightrefseq) = $ref3 =~ /^(.*$VRrefseq0)(.+)($VRrefseq1.*)$/; 
		my $VRpos30  = length($leftrefseq);
		my $VRrefseqlength  = length($VRrefseq3);
		my $VRpos31 = length($rightrefseq);
		($leftqueseq, $VRqueseq3, $rightqueseq) = $que3 =~ /^(.{${VRpos30}})(.+)(.{${VRpos31}})$/;

		$ref3Ind = 0;
		my ($VRpos20, $VRpos21) = (0,0);
		for (my $i = 0; $i < @ref2; $i++) {
			$ref3Ind ++ if $ref2 ne "-";
			if ($VRpos20 eq 0) {
				$VRpos20 = $i if $ref3Ind eq $VRpos30;
			}
			if ($VRpos21 eq 0) {
				$VRpos21 = $i if $ref3Ind eq $VRpos31;
			}
		}

		($leftqueseq2, $VRqueseq2, $rightqueseq2) = $que2 =~ /^(.{${VRpos20}})(.+)(.{${VRpos21}})$/;
		($leftrefseq2, $VRrefseq2, $rightrefseq2) = $ref2 =~ /^(.{${VRpos20}})(.+)(.{${VRpos21}})$/;


		my ($VRbadseq3, $VRins3, $VRdel3, $VRmat3, $VRmis3, $VRmetpos3, $VRmetneg3, $VRseqlength3) = get_seq_indel_stat($VRrefseq3, $VRqueseq3);
		my ($VRbadseq2, $VRins2, $VRdel2, $VRmat2, $VRmis2, $VRmetpos2, $VRmetneg2, $VRseqlength2) = get_seq_indel_stat($VRrefseq2, $VRqueseq2);

		my ($VRrefseqprint2, $VRqueseqprint2) = ($VRrefseq2, $VRqueseq2);
		my ($VRrefseqprint3, $VRqueseqprint3) = ($VRrefseq3, $VRqueseq3);

		my $VRmatperc2 = int(1000*$VRmat2 / $VRseqlength2+0.5)/10;
		my $VRmisperc2 = int(1000*$VRmis2 / $VRseqlength2+0.5)/10;
		my $VRinsperc2 = int(1000*$VRins2 / $VRseqlength2+0.5)/10;
		my $VRdelperc2 = int(1000*$VRdel2 / $VRseqlength2+0.5)/10;

		my $VRmatperc3 = int(1000*$VRmat3 / $VRseqlength3+0.5)/10;
		my $VRmisperc3 = int(1000*$VRmis3 / $VRseqlength3+0.5)/10;
		my $VRinsperc3 = int(1000*$VRins3 / $VRseqlength3+0.5)/10;
		my $VRdelperc3 = int(1000*$VRdel3 / $VRseqlength3+0.5)/10;
		
		$toprint .= "$read,GOOD,$linecount,all,VRpos2,$VRpos20,$VRpos21,$VRmetpos2,$VRmetneg2,$VRmatperc2,$VRmisperc2,$VRinsperc2,$VRdelperc2,$VRseqlength2\n";
		$toprint .= "$read,GOOD,$linecount,ref,VRseq2,$VRrefseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,que,VRseq2,$VRqueseqprint2\n";
		$toprint .= "$read,GOOD,$linecount,bad,VRseq2,$VRbadseq2\n";

		$toprint .= "$read,GOOD,$linecount,all,VRpos3,$VRpos30,$VRpos31,$VRmetpos3,$VRmetneg3,$VRmatperc3,$VRmisperc3,$VRinsperc3,$VRdelperc3,$VRseqlength3\n";
		$toprint .= "$read,GOOD,$linecount,ref,VRseq3,$VRrefseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,que,VRseq3,$VRqueseqprint3\n";
		$toprint .= "$read,GOOD,$linecount,bad,VRseq3,$VRbadseq3\n";

		@return = ($VRmetpos2, $VRmetneg2, $VRmatperc2, $VRmisperc2, $VRinsperc2, $VRdelperc2, $VRseqlength2, $toprint);

		my $ref2join = $VRrefseq2;
		my $que2join = $VRqueseq2;
		while ($ref2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $inslen = length($curr);
			my $pos21 = $pos20 + $inslen;
			$pos20 += $VRpos20;
			$pos21 += $VRpos20;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = "VR";
			$toprint .= "$read,GOOD,$linecount,big,ins,$regiontype,$pos20,$pos21,$ref2length,$inslen,$pos30,$pos31,$ref3length\n";
		}
		while ($que2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $dellen = length($curr);
			my $pos21 = $pos20 + $dellen;
			$pos20 += $VRpos20;
			$pos21 += $VRpos20;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = "VR";
			$toprint .= "$read,GOOD,$linecount,big,del,$regiontype,$pos20,$pos21,$ref2length,$dellen,$pos30,$pos31,$ref3length\n";
		}

		$ref2join = join("", @ref2);
		$que2join = join("", @que2);
		while ($ref2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $inslen = length($curr);
			my $pos21 = $pos20 + $inslen;
			$pos20 = $VRpos21 + 1 if ($pos20 <= $VRpos21 and $pos21 >  $VRpos21);
			$pos21 = $VRpos20 - 0 if ($pos20 <  $VRpos20 and $pos21 >= $VRpos20);
			next if $pos20 >= $VRpos20 and $pos20 <= $VRpos21 and $pos21 >= $VRpos20 and $pos21 <= $VRpos21;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = $pos30 < 100 ? "FWprimer" : $pos30 > @ref3 - 100 ? "RVprimer" : "NOTVR";
			$toprint .= "$read,GOOD,$linecount,big,ins,$regiontype,$pos20,$pos21,$ref2length,$inslen,$pos30,$pos31,$ref3length\n";
		}
		while ($que2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $dellen = length($curr);
			my $pos21 = $pos20 + $dellen;
			$pos20 = $VRpos21 + 1 if ($pos20 <= $VRpos21 and $pos21 >  $VRpos21);
			$pos21 = $VRpos20 - 0 if ($pos20 <  $VRpos20 and $pos21 >= $VRpos20);
			next if $pos20 >= $VRpos20 and $pos20 <= $VRpos21 and $pos21 >= $VRpos20 and $pos21 <= $VRpos21;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = $pos30 < 60 ? "FWprimer" : $pos30 > @ref3 - 60 ? "RVprimer" : "NOTVR";
			$toprint .= "$read,GOOD,$linecount,big,del,$regiontype,$pos20,$pos21,$ref2length,$dellen,$pos30,$pos31,$ref3length\n";
		}
	}

	else {
		my $ref2join = join("", @ref2);
		my $que2join = join("", @que2);
		while ($ref2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $inslen = length($curr);
			my $pos21 = $pos20 + $inslen;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = $pos30 < 60 ? "FWprimer" : $pos30 > @ref3 - 60 ? "RVprimer" : "NOTVR";
			$toprint .= "$read,GOOD,$linecount,big,ins,$regiontype,$pos20,$pos21,$ref2length,$inslen,$pos30,$pos31,$ref3length\n";
		}
		while ($que2join =~ /[\-]{4,999}/g) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my $pos20 = length($prev);
			my $dellen = length($curr);
			my $pos21 = $pos20 + $dellen;
			my $pos30 = $dict{ref2toref3}{$pos20};
			my $pos31 = $dict{ref2toref3}{$pos21};
			$pos30 = @ref3 if not defined $pos30;
			$pos31 = @ref3 if not defined $pos31;
			my $regiontype = $pos30 < 60 ? "FWprimer" : $pos30 > @ref3 - 60 ? "RVprimer" : "NOTVR";
			$toprint .= "$read,GOOD,$linecount,big,del,$regiontype,$pos20,$pos21,$ref2length,$dellen,$pos30,$pos31,$ref3length\n";
		}
		@return = (-1,-1,-1,-1,-1,-1,-1, $toprint);
	}
	return(@return);
}




sub get_seq_indel_stat {
	my ($refseq, $queseq) = @_;
	my @refseqs = split("", $refseq);
	my @queseqs = split("", $queseq);
	my ($badseq, $ins, $del, $mat, $mis, $metpos, $metneg) = ("", 0,0,0,0,0,0,0);
	my $length = 0;
	for (my $i = 0; $i < @refseqs; $i++) {
		my $refseq = $refseqs[$i];
		my $queseq = $queseqs[$i];
		if ($refseq ne "-" and $queseq ne "-" and $refseq eq $queseq) {
			$mat ++;
			$badseq .= ".";
			$length ++;
		}
		elsif ($refseq ne "-" and $queseq ne "-" and $refseq ne $queseq) {
			if ($refseq eq "G" and $queseq eq "A") {
				$metneg ++;
				$badseq .= ".";
			$length ++;
			}
			elsif ($refseq eq "C" and $queseq eq "T") {
				$metpos ++;
				$badseq .= ".";
				$length ++;
			}
			else {
				$mis ++;
				$badseq .= "m";
				$length ++;
			}
		}
		elsif ($refseq eq "N" and $queseq eq "-") {
			$badseq .= ".";
			$length ++;
		}
		elsif ($refseq eq "-" and $queseq ne "-") {
			$ins ++;
			$badseq .= "I";
			$length ++;
		}
		elsif ($refseq ne "-" and $queseq eq "-") {
			$del ++;
			$badseq .= "D";
			$length ++;
		}
		else {
			DIELOG($outLog, "Unexpected!\nrefseqs=" . join("", @refseqs) . "\nqueseqs=" . join("", @queseqs) . "\ni=$i, refseqnuc=$refseq, queseqnuc=$queseq\n\n");
		}
	}
	return($badseq, $ins, $del, $mat, $mis, $metpos, $metneg, $length);
}

sub check_file {
	my ($file, $type, $outLog) = @_;
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: $type file $file does not exist!\n") if ex($file) == 0;
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: $type file $file is empty!\n")       if -s $file  == 0;

	my $filetype = `file -b --mime-type $file`; chomp($filetype);
	my $cmd = ($file =~ /\.(rmdup|bam)$/ or $filetype =~ /(gzip|binary)/) ? "samtools view $file|" : "$file";
	my ($linecount, $check) = (0,0);
	my $currseq; 

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
		LOG($outLog, date() . "REF $def length=" . length($seq) . "\n","NA");
   }
   close $in;
	return(\%ref);
}

###########
# LOGFILE #
###########

sub parse_footLoop_logFile {
	my ($footLoop_folder, $outLog) = @_;
	$footLoop_folder = "footLoop_folder_unknown" if not defined $opt_n;


	# log file
	my $footLoop_logFile = "$footLoop_folder/logFile.txt";
	LOG($outLog, date() . "Logfile = $LRD$footLoop_logFile$N\n"); 

	# parse footLoop logfile
	DIELOG($outLog, "footLoop_2_filterBAMFile.pl: \n\nCan't find $footLoop_logFile! Please run footLoop.pl first before running this!\n\n") if not -e $footLoop_logFile;
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
		LOG($outLog, "\t$YW\$0${N}::${LGN}parse_bedFile$N: gene $gene = $length bp\n","NA");
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

