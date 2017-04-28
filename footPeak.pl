#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_g $opt_i $opt_p $opt_x $opt_y $opt_d $opt_s $opt_k);
getopts("vg:i:p:x:y:d:s:k:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/jeep/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;

my ($faFile, $indexFile, $peakFile, $x, $y, $min, $groupsize, $kmer) = ($opt_g, $opt_i, $opt_p, $opt_x, $opt_y, $opt_d, $opt_s, $opt_k);
die "usage: $0 -g hg19.fa -i UNMODIFIED geneindexes.bed -p CALM3_pos75.txt

Options (don't worry about this if you're just making length boxplots)
-x: left buffer (default: -10bp) 
-y: (right buffer) (default: +10bp)
-d: Maxmimum distance between 2 peaks (default: 150bp)
-s: Range of group size bin (default: 200bp)
-k: Kmer size +/- (default: 50bp)

" unless defined $faFile and defined $indexFile and defined $peakFile and -e $faFile and -e $indexFile and -e $peakFile;

$x = defined($x) ? $x : -10;
$y = defined($y) ? $y : 10;
$min = defined($min) ? $min : 150;
$groupsize = defined($groupsize) ? $groupsize : 200;
$kmer = defined($kmer) ? $kmer : 50;
print "X=$x, Y=$y\n";
#fastafrombed
system("bedtools_bed_change.pl -x $x -y $y -i $indexFile -o $indexFile\_$x\_$y.bed > /dev/null 2>&1");

# Get Real Coordinate of each Gene
my %coor;
my @LINE = `cat $indexFile\_$x\_$y.bed`;
print "\n$YW Parsing index file $indexFile\_$x\_$y.bed$N\n";
foreach my $line (@LINE) {
	chomp($line);
	my ($chr, $beg, $end, $name) = split("\t", $line);
	print "\tParsed:Gene=$CY$name$N, COORD=$CY$chr\t$beg\t$end\t$name$N\n\n";
	$coor{uc($name)}{chr} = $chr;
	$coor{uc($name)}{beg} = $beg;
	$coor{uc($name)}{end} = $end;
}

system("fastaFromBed -fi $faFile -bed $indexFile\_$x\_$y.bed -fo $indexFile\_$x\_$y.fa -name");
my $FA = $faFile;
$faFile = "$indexFile\_$x\_$y.fa";
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
#my ($folder1, $fileName1) = mitochy::getFilename($peakFile, "folder");

my %data; my $totalpeak = 0; my %end; my %bad;
open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
my ($gene) = $peakFile =~ /(\w+)_\w\w\w\d+(CG)?.txt/; 
die "Died cannot ge gene from peakfile $peakFile\n" unless defined $gene;
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
		if ($peak == 0 and $val[$i] == 9) {
			$peak = 1;
			$totalpeak ++;#= (keys %data);
			my $base = $seq{$gene}[$i]; die "Died i = $i name=$name total=${@val}\n" if not defined($base);
			my $index = int($i / $groupsize);
#			$data{$name}{$totalpeak}{beg} = $i;
#			push(@{$data{$name}{$totalpeak}{seq}}, $base);
			$data{$index}{$totalpeak}{beg} = $i;
			$data{$index}{$totalpeak}{name} = $name;
			$bad{$index}{$totalpeak}{beg} = 1 if $i <= 50;
			for (my $j = $i - 500; $j < $i; $j++) {
				next if $j < 0;
				last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
				$bad{$index}{$totalpeak}{beg} = 1 if $val[$j] == 6;
				last if $val[$j] == 6;
			}
			push(@{$data{$index}{$totalpeak}{seq}}, $base);
			$count = 0;
##			print "Peak $totalpeak: name $name: $base";
			$beg = $i; $end = 0;
			$peakcount ++;
		}
		elsif ($peak == 1 and (($val[$i] != 9 and $count == $min) or $i == $total - 1 or $i == @val - 1)) {
			my $base = $seq{$gene}[$i]; die "Died i = $i name=$name total=${@val}\n" if not defined($base);
			my $index = int($beg / $groupsize);
			$end = $i - $min;
			my $indexend = int($end/ $groupsize);
			push(@{$data{$index}{$totalpeak}{seq}}, $base);
			$bad{$index}{$totalpeak}{end} = 1 if $val[$end] == 6 or $i == $total - 1 or $i == @val - 1;
##			print "$base";
			$data{$index}{$totalpeak}{end} = $end;
			for (my $j = $end; $j < $j + 500; $j++) {
				last if $j == @val - 1 or $j == $total - 1;
				last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
				$bad{$index}{$totalpeak}{end} = 1 if $val[$j] == 6;
				last if $val[$j] == 6;
			}
			$end{$totalpeak} = $indexend;
##			print "\ni=$i indexBeg=$index indexEnd=$indexend Peak=$totalpeak name $name: BEG=$beg, END=$end\n\n";
#			push(@{$data{$name}{$totalpeak}{seq}}, $base);
#			$data{$name}{$totalpeak}{end} = $i;
			$count = 0;
			$peak = 0;
			$beg = 0; $end = 0;
		}
		elsif ($peak == 1 and $count < $min) {
			my $base = $seq{$gene}[$i]; 
			die "\nUndefined base at i=$i gene=$gene, total = $total\n\n" if not defined $base;
			$count ++ if $val[$i] != 9;
			$count = 0 if $val[$i] == 9;
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
			my $beg0 = $beg - $kmer < 0 ? 0 : $beg - $kmer;
			my $beg1 = $beg + $kmer > @seq-1 ? @seq-1 : $beg+$kmer;
			my $end0 = $end - $kmer < 0 ? 0 : $end - $kmer;
			my $end1 = $end + $kmer > @seq-1 ? @seq-1 : $end+$kmer;
			my $begSeq = join("", @seq[$beg0..$beg1]);
			my $endSeq = join("", @seq[$end0..$end1]);
			die "beg0=$beg0; beg1=$beg1; end0=$end0; end1=$end1\n" if not defined $seq[$beg0] or not defined($seq[$beg1]);
#			print "GENE=$gene\tGROUP=$indexBeg.$indexEnd\tBEG=$beg\tEND=$end\tPEAK=$peak\tNAME=$name\tbeg=$begSeq\tend=$endSeq\n";
		}
	}
}

open (my $outLen, ">", "$peakName.len") or die;
open (my $outOriginal, ">", "$peakName.original") or die;
my $strand = $peakFile =~ /Pos/i ? "pos" : $peakFile =~ /Neg/i ? "neg" : die "Cannot determine strand for file $peakFile\n";
open (my $out, ">", "$gene.prob");
my %kmer;
my $countz = 0;
foreach my $group (sort {$group{$b} <=> $group{$a} || $a cmp $b} keys %group) {
	$countz ++;
#	printf "Group=$countz (name=$group)\tTotalPeak=$group{$group}\tPercent=%.1f\n", $group{$group}/$totalPeak * 100;
#	next if $group{$group} / $totalPeak * 100 < 5;
	my ($indexBeg, $indexEnd) = $group =~ /^(\d+)\.(\d+)$/;
	foreach my $peak (sort {$a <=> $b} keys %{$final{$indexBeg}{$indexEnd}}) {
		my $seq  = $final{$indexBeg}{$indexEnd}{$peak}{seq};
		my $name = $final{$indexBeg}{$indexEnd}{$peak}{name};
		my $beg  = $final{$indexBeg}{$indexEnd}{$peak}{beg};
		my $end  = $final{$indexBeg}{$indexEnd}{$peak}{end};
		my ($gene) = $name =~ /^(.+)\.SEQ_/;
		my ($origChr, $origBeg, $origEnd) = ($coor{$gene}{chr}, $coor{$gene}{beg}, $coor{$gene}{end});
		my $peakBeg = $origBeg + $beg - 1; $peakBeg = $beg if $peakBeg < $beg;
		my $peakEnd = $origBeg + $end - 1; $peakEnd = $beg if $peakEnd < $beg;
		#die "Gene = $gene, name = $name, CHR=$origChr, BEG=$origBeg, END=$origEnd, beg = $beg, end = $end, peakBeg = $peakBeg, peakEnd = $peakEnd\n";
		my $len = $end - $beg + 1;
		my $newstrand = $strand eq "pos" ? "+" : $strand eq "neg" ? "-" : ".";
		print $outLen "$gene\t$beg\t$end\t$name\t$len\t$newstrand\n";
		print $outOriginal "$origChr\t$peakBeg\t$peakEnd\t$name\t$len\t$newstrand\n";
		my $badbeg  = $final{$indexBeg}{$indexEnd}{$peak}{badbeg};
		my $badend  = $final{$indexBeg}{$indexEnd}{$peak}{badend};
		my $beg0 = $beg - $kmer < 0 ? 0 : $beg - $kmer;
		my $beg1 = $beg + $kmer > @seq-1 ? @seq-1 : $beg+$kmer;
		my $end0 = $end - $kmer < 0 ? 0 : $end - $kmer;
		my $end1 = $end + $kmer > @seq-1 ? @seq-1 : $end+$kmer;
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
		print "\n\n!!!!!!!!!!!!!!!!\nCONTACT STELLA THIS HAPPEN!!! AT: GENE $gene GROUP $group PEAK $peak SEQ $seq NAME $name beg0=$beg0; beg1=$beg1; end0=$end0; end1=$end1\n!!!!!!!!!!!!!!!!!!!!\n" if not defined $seq[$beg0] or not defined($seq[$beg1]);
#		die "GROUP=$countz (name=$group)\tGENE=$gene\tGROUP=$indexBeg.$indexEnd\tBEG=$beg\tEND=$end\tPEAK=$peak\tNAME=$name\tbeg=$begSeq\tend=$endSeq\n";
	}
}
close $outLen;

print "\nOUTPUT:\n\t- Relative coordinate bed file: $LGN$peakName.len$N\n";
print "\t- UCSC Original coordinate bed file: $LGN$peakName.original$N\n\n--------------------\nKMER:";
#open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
#close $out1;
my @seqs = qw(A C G T);
print "type";
foreach my $seq1(sort @seqs) {
	foreach my $seq2 (sort @seqs) {
		print "\t$seq1$seq2";
	}
}
print "\n";
kmer(join("", @seq),"all");

my @type = qw(beg mid end);
my $obs; my $exp; my $odd;
foreach my $type (@type[0..2]) {
	$obs .= "$type";
	$exp .= "any" if $type eq "beg";
	$odd .= "$type";
	foreach my $seq1(sort @seqs) {
		foreach my $seq2 (sort @seqs) {
			my $kmer = "$seq1$seq2";
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
#			print "\tcount=$count;total=$total;perc=$perc";
			$obs .= "\t$perc";
#			$obs .= "\tcount=$count;total=$total;perc=$perc";
			$exp .= "\t$percE" if $type eq "beg";
#			$exp .= "\tcount=$countE;total=$totalE;perc=$percE" if $type eq "beg";
			$odd .= "\t$ratio";
		}
	}
	$obs .= "\n";
	$exp .= "\n" if $type eq "beg";
	$odd .= "\n";
}

print ">Observed:\n$obs\n\n";
print ">Expected:\n$exp\n\n";
print ">OddRatio:\n$odd\n\n";
sub kmer {
	my ($seq, $type) = @_;
	$seq = uc($seq);
	my @seqz = split("", $seq);
	for (my $i = 0; $i < @seqz-1; $i++) {
		my $kmer = "$seqz[$i]$seqz[$i+1]";
		$kmer{$type}{$kmer} ++;
		$kmer{$type}{total} ++;
	}
#	for (my $i = 0; $i < @seqz-1; $i++) {
#		my $kmer = "$seqz[$i]$seqz[$i+1]";
#		print "$type $kmer = $kmer{$type}{$kmer}\n";
#	}
}
