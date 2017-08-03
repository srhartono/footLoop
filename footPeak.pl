#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_g $opt_i $opt_p $opt_x $opt_y $opt_d $opt_s $opt_k $opt_K $opt_n $opt_h);
getopts("vg:i:p:x:y:d:s:k:K:n:h");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $die = "\nDied at file$CY" . __FILE__ . "$N at line $LGN";
my ($faFile, $indexFile, $peakFile, $x, $y, $min, $groupsize, $dist, $outDir) = ($opt_g, $opt_i, $opt_p, $opt_x, $opt_y, $opt_d, $opt_s, $opt_k, $opt_n);
my $usage = "
Usage: $YW$0$N -g$CY <genomic fasta>$N -i$LPR <UNMODIFIED geneIndexes.bed>$N -p$LGN <Peak file>$N -n$YW <Output Folder>$N

${YW}Example: $YW$0$N -g$CY hg19.fa$N -i$LPR geneIndexes.bed$N -p$LGN CALM3_Pos75.txt$N -n$YW CALM3_Pos75_Out$N

$LRD	========== !!IMPORTANT!! ========== $N
1.	The format of the file from -p *has* to be: ${CY}GENE$N\_NNNDD.txt
	NNN is Pos or Neg
	DD is the percent threshold (e.g. 75)
	If there's 'CG' after DD that's okay!
	E.g.:$YW CALM3_pos75.txt or CALM3_NEG75.txt or CALM3_Pos75CG.txt$N

2.	The GENE will  be used to extract gene from geneindexes.bed (case insensitive!)
	So if in geneindexes.bed, the gene name is CALM and the file name is CALM3_pos40.txt,
	this will >not< work as CALM is >not<the same as CALM3


3.	Put -x and -y EXACTLY like what you used for footLoop.pl
	e.g. in footLoop.pl you gave -100 left buffer and 100 right buffer
	do this: -x -10 -y 10
	If you didn't specify -x and -y when doing footLoop.pl, then don't put anything here as well

	${LRD}The result will be wrong if -x and -y isn't the same as what you used in footLoop.pl!$N
";

my $usage_long = "
$LRD	=================================== $N

${LGN}Options$N [default] 
-x: Length of left 'buffer' in basepair [0]
-y: Length of right 'buffer' in basepair [0]
-d: Maxmimum distance between 2 peaks in basepair [150]
-s: Range of group size bin in basepair [200]
-k: Length of seequence to take from start/end of each peak [50]
-K: [2] Kmer size; the 'k' of 'k'mer [-K 2 (AA/AT/AG/etc];
	 -K 2 -> AA/AC/AG/AT/CA/CC/...
	 -K 3 -> AAA/AAC/AAG/AAT/ACA/ACT/...
	 -K 4 -> AAAA/AAAC/AAAG/AAAT/...

${LGN}Determining cluster$N

Example:

...|100      |110      |120 ... |250
...01234567${LGN}8${N}9012345678901234...901234${LGN}5${N}6789012345678....   <=$LGN this is position 1 at 0 means 100, 1 at 250 means 251$N
-----------|         PEAK   ...      |-----------

Above, peak begin at position 108 and end at position 255.

-d of 150:
Say there's another peak >> WITHIN SAME READ << that begins at position 280 and ending at position 500.
This peak is 'too close' to the first peak as its beginning (280) is less than 150bp away from the first peak's end (255)
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

";

die $usage . $usage_long if defined $opt_h;
die $usage unless defined $faFile and defined $indexFile and defined $peakFile and -e $faFile and -e $indexFile and -e $peakFile and defined $outDir;

my ($K) = defined ($opt_K) ? $opt_K : 2;
die $die . __LINE__ . "$N: -K *must* be 2 or 3 or 4! (Currently:$LGN$K$N)\n\n" unless $K =~ /^[234]$/;
my @kmer = @{make_kmer($K)};

if (not -d $outDir) {
	mkdir $outDir or die $die . __LINE__ . "$N: Cannot create directory (-n) $outDir: $!\n\n";
}

$x = defined($x) ? $x : 0;
$y = defined($y) ? $y : 0;
$min = defined($min) ? $min : 150;
$groupsize = defined($groupsize) ? $groupsize : 200;
$dist = defined($dist) ? $dist : 50;

# Get Real Coordinate of each Gene
my %coor;
$faFile = parse_index($indexFile, $x, $y);

# get seq from faFile
open (my $in, "<", $faFile) or die;
my %seq;
my $def; my @seq;
while (my $line = <$in>) {
	chomp($line);
	if ($line =~ /^>/) {
		if (defined($def)) {
			@{$seq{$def}} = @seq;
			@seq = (); undef $def;
		}
		($def) = $line =~ /^>(.+)$/;
		$def = uc($def);
		#print "Def = $def\n";
	}
	else {
		@seq = @seq == 0 ? split("", $line) : (@seq, split("", $line));
	}
}
@{$seq{$def}} = @seq;
close $in;

# process peak

#PEAK	.......999..9.999.....{300 more dots}...999.999999........
#           | P1| | P2 |                      | Peak 3  |         # if max non-9 buffer length is < 2
#           |  Peak1   |                      | Peak 2  |         # if max non-9 buffer length is < 300
#           |              Peak 1                       |         # if max non-9 buffer length is < 500

my %data; my $totalpeak = 0; my %end; my %bad;
open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
my ($gene) = $peakFile =~ /(\w+)_\w\w\w\d+(CG)?.txt/; 
die "Died cannot get gene from peakfile $peakFile\n" unless defined $gene;
$gene = uc($gene);
my ($peakName) = $peakFile =~ /(\w+_\w\w\w\d+(CG)?).txt/; $gene = uc($gene);
my $total = @{$seq{$gene}}; my $linecount = 0; my $peakcount = 0;
#print "Gene = $gene total = $total\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my ($name, @val) = split("\t", $line);
	$name = "$gene.$name";
	print "Line error skipped: $line\n" and next if not defined($name) or @val == 0;
#	print "\n\n>$name\n";	
	my $peak = 0; my $count = 0; my $beg = 0; my $end = 0;
	for (my $i = 0; $i < @val; $i++) {
#		next if $val[$i] eq "";
	#	print "DIED\n\n$line\n\nVal:" . join(" ", @val[@val-10..@val-1]) . "\n" and die if $val[$i] eq "";
		# 1. Not currently at peak, and current is peak (9)
		if ($peak == 0 and $val[$i] == 9) {
			$peak = 1;
			$totalpeak ++;#= (keys %data);
			my $base = $seq{$gene}[$i]; die "Died i = $i name=$name total=${@val}\n" if not defined($base);
			my $index = int($i / $groupsize);
			$data{$index}{$totalpeak}{beg} = $i;
			$data{$index}{$totalpeak}{name} = $name;
			$bad{$index}{$totalpeak}{beg} = 1 if $i <= 50;
			my $j0 = $i < 500 ? 0 : $i - 500;
			for (my $j = $j0; $j < $i; $j++) {
				last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
				$bad{$index}{$totalpeak}{beg} = 1 if $val[$j] == 6;
				last if $val[$j] == 6;
			}
			push(@{$data{$index}{$totalpeak}{seq}}, $base);
			$count = 0;
##			print "Peak $totalpeak: name $name: $base";
			$beg = $i; $end = $i;
			$peakcount ++;
		}
		# 2. Currently is peak (peak == 1) but value isn't 9 and total non-peak is equal to maximum buffer length, then peak is end
		elsif ($peak == 1 and (($val[$i] != 9 and $count == $min) or $i == $total - 1 or $i == @val - 1)) {
			my $base = $seq{$gene}[$i]; die "Died i = $i name=$name total=${@val}\n" if not defined($base);
			my $index = int($beg / $groupsize);
			#$end = $i - $min;
			my $indexend = int($end/ $groupsize);
			push(@{$data{$index}{$totalpeak}{seq}}, $base);
			$bad{$index}{$totalpeak}{end} = 1 if $val[$end] == 6 or $i == $total - 1 or $i == @val - 1;
			$data{$index}{$totalpeak}{end} = $end;
#			for (my $j = $end; $j < $j + 500; $j++) {
			for (my $j = $end; $j < $end + 500; $j++) {
				last if $j == @val - 1 or $j == $total - 1;
				last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
				$bad{$index}{$totalpeak}{end} = 1 if $val[$j] == 6;
				last if $val[$j] == 6;
			}
			$end{$totalpeak} = $indexend;
##			print "\ni=$i indexBeg=$index indexEnd=$indexend Peak=$totalpeak name $name: BEG=$beg, END=$end\n\n";
			$count = 0;
			$peak = 0;
			$beg = 0; $end = 0;
		}
		# 3. Currently is peak (peak == 1), but total non-peak length is less than required maximum non-peak length ($min)
		elsif ($peak == 1 and $count < $min) {
			my $base = $seq{$gene}[$i]; 
			die "\nUndefined base at i=$i gene=$gene, total = $total\n\n" if not defined $base;
			$count ++ if $val[$i] != 9;
			$count = 0 if $val[$i] == 9;
			$end = $i if $val[$i] == 9;
			my $index = int($beg / $groupsize);
			push(@{$data{$index}{$totalpeak}{seq}}, $base);
##			print "$base";
		}
	#	print "\nBEG=$beg, END=$end\n" if $i == @val - 1 or $i == $total - 1;
		last if $i == $total - 1;
	}
	$linecount ++;
##	last if $totalpeak > 2;
}
close $in1;
print "Processed $linecount lines, peak = $peakcount\n";
my %final;
##print "\n\n";
foreach my $index (sort {$a <=> $b} keys %data) {
	foreach my $peaknum (sort {$a <=> $b} keys %{$data{$index}}) {
		my $seq = join("", @{$data{$index}{$peaknum}{seq}});
		my $name = $data{$index}{$peaknum}{name};
		my $indexend = $end{$peaknum};
		my $badbeg = defined($bad{$index}{$totalpeak}{beg}) ? $bad{$index}{$totalpeak}{beg} : 0;
		my $badend = defined($bad{$index}{$totalpeak}{end}) ? $bad{$index}{$totalpeak}{end} : 0;
		$final{$index}{$indexend}{$peaknum}{badbeg} = $badbeg;
		$final{$index}{$indexend}{$peaknum}{badend} = $badend;
		$final{$index}{$indexend}{$peaknum}{seq} = $seq;
		$final{$index}{$indexend}{$peaknum}{name} = $name;
		$final{$index}{$indexend}{$peaknum}{beg} = $data{$index}{$peaknum}{beg};
		$final{$index}{$indexend}{$peaknum}{end} = $data{$index}{$peaknum}{end};
		
##		print "$index $indexend Peak=$peaknum $name $seq\n";
	}
}
@seq = @{$seq{$gene}};
my %group; my $totalPeak = 0;
foreach my $indexBeg (sort {$a <=> $b} keys %final) {
	foreach my $indexEnd (sort {$a <=> $b} keys %{$final{$indexBeg}}) {
		foreach my $peak (sort {$a <=> $b} keys %{$final{$indexBeg}{$indexEnd}}) {
			$group{"$indexBeg.$indexEnd"} ++;
			$totalPeak ++;
			my $seq  = $final{$indexBeg}{$indexEnd}{$peak}{seq};
			my $name = $final{$indexBeg}{$indexEnd}{$peak}{name};
			my $beg  = $final{$indexBeg}{$indexEnd}{$peak}{beg};
			my $end  = $final{$indexBeg}{$indexEnd}{$peak}{end};
			my $badbeg  = $final{$indexBeg}{$indexEnd}{$peak}{badbeg};
			my $badend  = $final{$indexBeg}{$indexEnd}{$peak}{badend};
			my $beg0 = $beg - $dist < 0 ? 0 : $beg - $dist;
			my $beg1 = $beg + $dist > @seq-1 ? @seq-1 : $beg+$dist;
			my $end0 = $end - $dist < 0 ? 0 : $end - $dist;
			my $end1 = $end + $dist > @seq-1 ? @seq-1 : $end+$dist;
			my $begSeq = join("", @seq[$beg0..$beg1]);
			my $endSeq = join("", @seq[$end0..$end1]);
			die "beg0=$beg0; beg1=$beg1; end0=$end0; end1=$end1\n" if not defined $seq[$beg0] or not defined($seq[$beg1]);
#			print "GENE=$gene\tGROUP=$indexBeg.$indexEnd\tBEG=$beg\tEND=$end\tPEAK=$peak\tNAME=$name\tbeg=$begSeq\tend=$endSeq\n";
		}
	}
}

open (my $outLen, ">", "$outDir/$peakName.len") or die "Cannot write to $outDir/$peakName.len: $!\n";
open (my $outOriginal, ">", "$outDir/$peakName.original") or die "Cannot write to $outDir/$peakName.original: $!\n";
my $strand = $peakFile =~ /Pos/i ? "pos" : $peakFile =~ /Neg/i ? "neg" : die "Cannot determine strand for file $peakFile\n";
open (my $outz, ">", "$outDir/$gene.group") or die "Cannot write to $outDir/$gene.group: $!\n";
my %kmer;
my $countz = 0;
foreach my $group (sort {$group{$b} <=> $group{$a} || $a cmp $b} keys %group) {
	$countz ++;
#	next if $group{$group} / $totalPeak * 100 < 5;
	printf $outz "Group=$countz (name=$group)\tTotalPeak=$group{$group}\tPercent=%.1f\n", $group{$group}/$totalPeak * 100;
	my ($begPrint, $midPrint, $endPrint) = ("","","");
	my ($indexBeg, $indexEnd) = $group =~ /^(\d+)\.(\d+)$/;
	foreach my $peak (sort {$a <=> $b} keys %{$final{$indexBeg}{$indexEnd}}) {
		my $seq  = $final{$indexBeg}{$indexEnd}{$peak}{seq};
		my $name = $final{$indexBeg}{$indexEnd}{$peak}{name};
		my $beg  = $final{$indexBeg}{$indexEnd}{$peak}{beg};
		my $end  = $final{$indexBeg}{$indexEnd}{$peak}{end};
		my ($gene) = $name =~ /^(.+)\.SEQ_/;
		my ($origChr, $origBeg, $origEnd) = ($coor{$gene}{chr}, $coor{$gene}{beg}, $coor{$gene}{end});
		my $peakBeg = $origBeg + $beg; $peakBeg = $origBeg if $peakBeg < $origBeg;
		my $peakEnd = $origBeg + $end + 1; die "$peak: start ($peakBeg) is less than end ($peakEnd)\n" if $peakEnd < $peakBeg;
		#die "Gene = $gene, name = $name, CHR=$origChr, BEG=$origBeg, END=$origEnd, beg = $beg, end = $end, peakBeg = $peakBeg, peakEnd = $peakEnd\n";
		my $len = $end - $beg + 1;
		my $newstrand = $strand eq "pos" ? "+" : $strand eq "neg" ? "-" : ".";
		print $outLen "$gene\t$beg\t$end\t$name\t$len\t$newstrand\n";
		print $outOriginal "$origChr\t$peakBeg\t$peakEnd\t$name\t$len\t$newstrand\n";
		my $badbeg  = $final{$indexBeg}{$indexEnd}{$peak}{badbeg};
		my $badend  = $final{$indexBeg}{$indexEnd}{$peak}{badend};
		my $beg0 = $beg - $dist < 0 ? 0 : $beg - $dist;
		my $beg1 = $beg + $dist > @seq-1 ? @seq-1 : $beg+$dist;
		my $end0 = $end - $dist < 0 ? 0 : $end - $dist;
		my $end1 = $end + $dist > @seq-1 ? @seq-1 : $end+$dist;
		my $begSeq = $badbeg == 0 ? join("", @seq[$beg0..$beg1]) : "NA";
		my $endSeq = $badend == 0 ? join("", @seq[$end0..$end1]) : "NA";
		my $temp = $begSeq;
		my $midSeq = join("", @seq[$beg1..$end0]);
		$begSeq = uc(reverse($endSeq)) if $strand eq "neg";
		$begSeq =~ tr/ATGC/TACG/ if $strand eq "neg";
		$endSeq = uc(reverse($temp)) if $strand eq "neg";
		$endSeq =~ tr/ATGC/TACG/ if $strand eq "neg";
		$midSeq = uc(reverse($midSeq)) if $strand eq "neg";
		$midSeq =~ tr/ATGC/TACG/ if $strand eq "neg";
		kmer($begSeq, "beg");
		kmer($endSeq, "end");
		kmer($midSeq, "mid");
		$begPrint .= "$gene\t$beg0\t$beg1\t$peak.beg\t0\t$newstrand\t$begSeq\n";
		$midPrint .= "$gene\t$beg1\t$end0\t$peak.mid\t0\t$newstrand\t$midSeq\n";
		$endPrint .= "$gene\t$end0\t$end1\t$peak.end\t0\t$newstrand\t$endSeq\n";
		print "\n\n!!!!!!!!!!!!!!!!\nCONTACT STELLA THIS HAPPEN!!! AT: GENE $gene GROUP $group PEAK $peak SEQ $seq NAME $name beg0=$beg0; beg1=$beg1; end0=$end0; end1=$end1\n!!!!!!!!!!!!!!!!!!!!\n" if not defined $seq[$beg0] or not defined($seq[$beg1]);
#		die "GROUP=$countz (name=$group)\tGENE=$gene\tGROUP=$indexBeg.$indexEnd\tBEG=$beg\tEND=$end\tPEAK=$peak\tNAME=$name\tbeg=$begSeq\tend=$endSeq\n";
	}
	print $outz "\n$begPrint$midPrint$endPrint";
}
close $outLen;
close $outz;

print "\nOUTPUT:\n";
print "\t- Relative coordinate bed file: $LGN$outDir/$peakName.len$N\n";
print "\t- UCSC Original coordinate bed file: $LGN$outDir/$peakName.original$N\n";
print "\t- Kmer odds ratio file: $LGN$outDir/$peakName.kmer$N\n";
#print "\n\n--------------------\nKMER:";
#print "\n\n--------------------\n";

open (my $outKmer, ">", "$outDir/$peakName.kmer") or die "Cannot write to $outDir/$peakName.kmer: $!\n";
my @seqs = @{make_kmer($K)};
print $outKmer "type";
foreach my $seq1(sort @seqs) {
	print $outKmer "\t$seq1";
}
print $outKmer "\n";
kmer(join("", @seq),"all");

my @type = qw(beg mid end);
my $obs; my $exp; my $odd;
foreach my $type (@type[0..2]) {
	$obs .= "$type";
	$exp .= "any" if $type eq "beg";
	$odd .= "$type";
	foreach my $seq1(sort @seqs) {
		my $kmer = "$seq1";
		my $count = defined $kmer{$type}{$kmer} ? $kmer{$type}{$kmer} : 0;
		my $total = defined $kmer{$type}{total} ? $kmer{$type}{total} : 0;
		my $countE = defined $kmer{all}{$kmer} ? $kmer{all}{$kmer} : 0;
		my $totalE = defined $kmer{all}{total} ? $kmer{all}{total} : 1;
		my $perc  = int((1+$count*1000)/(1+$total))/10;
		my $percE = int((1+$countE*1000)/(1+$totalE))/10;
		my $ratio = $percE == 0 ? 0 : int(100*$perc / $percE)/100;
		$perc = $perc <= 2 ? "$LBU$perc$N" : $perc <= 5 ? "$LCY$perc$N" : $perc <= 7 ? "$N$perc$N" : $perc <= 10 ? "$YW$perc$N" : "$LRD$perc$N";
		$percE = $percE <= 2 ? "$LBU$percE$N" : $percE <= 5 ? "$LCY$percE$N" : $percE <= 7 ? "$N$percE$N" : $percE <= 10 ? "$YW$percE$N" : "$LRD$percE$N";
		$ratio = $ratio <= 0.5 ? "$LBU$ratio$N" : $ratio <= 0.9 ? "$LCY$ratio$N" : $ratio <= 1/0.9 ? "$N$ratio$N" : $ratio <= 1/0.5 ? "$YW$ratio$N" : "$LRD$ratio$N";
#		print $outKmer "\tcount=$count;total=$total;perc=$perc";
		$obs .= "\t$perc";
#		$obs .= "\tcount=$count;total=$total;perc=$perc";
		$exp .= "\t$percE" if $type eq "beg";
#		$exp .= "\tcount=$countE;total=$totalE;perc=$percE" if $type eq "beg";
		$odd .= "\t$ratio";
	}
	$obs .= "\n";
	$exp .= "\n" if $type eq "beg";
	$odd .= "\n";
}

print $outKmer ">Observed:\n$obs\n\n";
print $outKmer ">Expected:\n$exp\n\n";
print $outKmer ">OddRatio:\n$odd\n\n";
sub kmer {
	my ($seq, $type) = @_;
	$seq = uc($seq);
	my @seqz = split("", $seq);
	for (my $i = 0; $i < @seqz-($K-1); $i++) {
		my $kmer;
		for (my $j = $i; $j < $i+$K; $j++) {
			$kmer .= $seqz[$j];
		}
		$kmer{$type}{$kmer} ++;
		$kmer{$type}{total} ++;
	}
#	for (my $i = 0; $i < @seqz-1; $i++) {
#		my $kmer = "$seqz[$i]$seqz[$i+1]";
#		print "$type $kmer = $kmer{$type}{$kmer}\n";
#	}
}

sub parse_index {
	my ($indexFile, $leftBuf, $riteBuf) = @_;
	system("bedtools_bed_change.pl -x $leftBuf -y $riteBuf -i $indexFile -o $indexFile\_$leftBuf\_$riteBuf.bed > /dev/null 2>&1");

	my @LINE = `cat $indexFile\_$leftBuf\_$riteBuf.bed`;
	print "\n$YW Parsing index file $indexFile\_$leftBuf\_$riteBuf.bed$N\n";
	foreach my $line (@LINE) {
		chomp($line);
		my ($chr, $beg, $end, $name) = split("\t", $line);
		print "\tParsed:Gene=$CY$name$N, COORD=$CY$chr\t$beg\t$end\t$name$N\n\n";
		$coor{uc($name)}{chr} = $chr;
		$coor{uc($name)}{beg} = $beg;
		$coor{uc($name)}{end} = $end;
	}
	system("fastaFromBed -fi $faFile -bed $indexFile\_$leftBuf\_$riteBuf.bed -fo $indexFile\_$leftBuf\_$riteBuf.fa -name");
	$faFile = "$indexFile\_$leftBuf\_$riteBuf.fa";
	return($faFile);
}

sub make_kmer {
	# Initialize Count tables. Don't mind the commented-out, these are for debug.
	my ($order) = @_;
	my @preword = qw(A C G T);
	my @words = qw(A C G T);
	for (my $i = 1; $i < $order; $i++) {
		my @curr;
		for (my $k = 0; $k < @preword; $k++) {
			for (my $j = 0; $j < @words; $j++) {
				push(@curr, "$preword[$k]$words[$j]");
			}
		}
		@preword = @curr;
	}
	my $max = @preword > 9 ? 9 : @preword-1;
	#print "\nK=$K; Kmer = " . join(",", @preword[0..$max]);
	#print "\n" if @preword <= 9;
	#print "..." . "(total =$LGN " . scalar(@preword) . "$N)\n\n" if @preword > 9;
	return(\@preword);
}

__END__
	my %count;
	for (my $i = 0; $i < @preword; $i++) {
		for (my $j = 0; $j < @words; $j++) {
			$count{$preword[$i]}{$words[$j]} = 0;
		}
	}
#	print "\n";
	return %count;
}