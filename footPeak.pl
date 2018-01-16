#!/usr/bin/perl

# 0 is bad or not data (6)
# 2 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak

use strict; use warnings; use Getopt::Std; use Time::HiRes; use Benchmark qw(:all); use Benchmark ':hireswallclock'; use Carp;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_g $opt_i $opt_p $opt_x $opt_y $opt_d $opt_s $opt_k $opt_K $opt_n $opt_h $opt_t $opt_w $opt_l);
getopts("vg:i:p:x:y:d:s:k:K:n:ht:w:l:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';

my ($footFolder) = $opt_i;
my $date = getDate();
my $uuid = getuuid();

my ($usage) = check_sanity($footFolder);
my $logFile = "$footFolder/logFile.txt";
my ($opts, $outLog) = parse_footLoop_logFile($logFile, $date, $uuid);

my $faFile        = $opts->{g};
my $indexFile		= $opts->{i};
my $origFolder 	= $opts->{origDir};
my $x 				= $opts->{x};
my $y 				= $opts->{y};
my $threshold 		= $opts->{t};
my $window 			= $opts->{w};
my $min 				= $opts->{d};
my $groupsize 		= $opts->{s};
my $minDis 			= $opts->{k};
my $outDir 			= $opts->{n};
my $Kmerz			= $opts->{K};
print $outLog . "$N: -K *must* be 2 or 3 or 4! (Currently:$LGN$Kmerz$N)\n\n" and die unless $Kmerz =~ /^[234]$/;

# Create kmer (deprecated)
my @kmer = @{make_kmer($Kmerz)};

#my ($faFile, $indexFile, $x, $y, $min, $groupsize, $minDis, $outDir) = ($opts->{g}, $opts->{i}, $opts->{x}, $opts->{y}, $opt_d, $opt_s, $opt_k, $opt_n);
$outDir = $origFolder if not defined $outDir;
print $outLog $usage . "faFile=$faFile\nindexFile=$indexFile\norigFolder=$origFolder\noutdir=$outDir\n" and die unless defined $faFile and defined $indexFile and defined $origFolder and -e $faFile and -e $indexFile and -e $origFolder and defined $outDir;
#die $usage if mydefined($faFile, $indexFile, $x, $y, $min, $groupsize, $minDis, $outDir) != 0;
if (not -d $outDir) {
	mkdir $outDir or die DIE() . __LINE__ . "$N: Cannot create directory (-n) $outDir: $!\n\n";
}

# Get Real Coordinate of each Gene
my %coor;
$indexFile = parse_index($indexFile, $outDir, $x, $y);
my @indexLine = `cat $indexFile`;
my $SEQ;
foreach my $line (@indexLine) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	$gene = uc($gene);
	$SEQ->{$gene} = "$chr\t$beg\t$end\t$gene\t$zero\t$strand";
	print "gene=$gene, coor=$chr, $beg, $end, $gene, $zero, $strand\n";
}
# get seq from faFile
open (my $in, "<", $faFile) or die;
my $fasta = new FAlite($in);
my %seq;
my $def; my @seq;
while (my $entry = $fasta->nextEntry()) {
	my $def = $entry->def; $def =~ s/^>//;
	my @seq = split("", $entry->seq);
	$seq{$def}{seq} = \@seq;
	$seq{$def}{loc} = findCGPos(\@seq);
	next;
	my ($C, $G, $tot) = (0,0,0);
	for (my $i = 0; $i < @seq-100; $i+= 10) {
		my $C1 = join("", @seq[$i..$i+100]) =~ tr/Cc/Cc/;
		my $G1 = join("", @seq[$i..$i+100]) =~ tr/Gg/Gg/;
		$C += $C1; $G += $G1;
		$tot ++;
	}
	$C = int($C/$tot*10)/10;
	$G = int($G/$tot*10)/10;
	print "$def: C=$C, G=$G\n";
}
close $in;
# process peak

#PEAK	.......999..9.999.....{300 more dots}...999.999999........
#           | P1| | P2 |                      | Peak 3  |         # if max non-9 buffer length is < 2
#           |  Peak1   |                      | Peak 2  |         # if max non-9 buffer length is < 300
#           |              Peak 1                       |         # if max non-9 buffer length is < 500

my @origFile = <$origFolder/*.orig>;

for (my $i = 0; $i < @origFile; $i++) {
	my $peakFile = $origFile[$i];
#debug
	next if $peakFile !~ /CALM3_Pos/;
##
	my ($peakFolder, $peakFilename) = getFilename($peakFile, "folder");
	$peakFilename =~ s/.orig$//;
	my ($gene, $strand) = $peakFilename =~ /^(\w+)_(Pos|Neg|Unk)$/; $gene = uc($gene);
	my ($totalPeak, $linecount, $peakcount, $total, %data, %end, %bad, %final, %group) = (0,0,0, scalar(@{$seq{$gene}{seq}}));
	print "FILENAME=$peakFilename, GENE=$gene\n";
	print $outLog "Died cannot get gene from peakfile $peakFile\n" and die unless defined $gene;
	print $outLog "Cannot find sequence of gene=$LCY$gene$N in $faFile!\n" and die if not defined $seq{$gene};
	my %pk = ("CH" => 0, "CG" => 0, "GH" => 0, "GC" => 0);
	print "$LGN$i$N. peakFile=$LCY$peakFile$N, Gene = $LRD$gene$N total length = $LPR$total$N\n";
	open (my $outPEAKCG, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_CG.PEAK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_CG.PEAK: $!\n";
	open (my $outNOPKCG, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_CG.NOPK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_CG.NOPK: $!\n";
	open (my $outPEAKCH, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_CH.PEAK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_CH.PEAK: $!\n";
	open (my $outNOPKCH, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_CH.NOPK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_CH.NOPK: $!\n";
	open (my $outPEAKGC, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_GC.PEAK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_GC.PEAK: $!\n";
	open (my $outNOPKGC, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_GC.NOPK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_GC.NOPK: $!\n";
	open (my $outPEAKGH, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_GH.PEAK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_GH.PEAK: $!\n";
	open (my $outNOPKGH, ">", "$peakFolder/$peakFilename\_$window\_$threshold\_GH.NOPK") or die "Failed to write to $peakFolder/$peakFilename\_$window\_$threshold\_GH.NOPK: $!\n";
	open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
	my ($l0, $t0) = (0,Benchmark->new());
	print $outPEAKCH "$peakFile\tPEAK\t$gene\tCH\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outNOPKCH "$peakFile\tNOPK\t$gene\tCH\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outPEAKCG "$peakFile\tPEAK\t$gene\tCG\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outNOPKCG "$peakFile\tNOPK\t$gene\tCG\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outPEAKGH "$peakFile\tPEAK\t$gene\tGH\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outNOPKGH "$peakFile\tNOPK\t$gene\tGH\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outPEAKGC "$peakFile\tPEAK\t$gene\tGC\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	print $outNOPKGC "$peakFile\tNOPK\t$gene\tGC\t$strand\t" . join("\t", @{$seq{$gene}{seq}}) . "\n";
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		next if $line =~ /^#/;
		my ($name, $type, $val) = split("\t", $line);
		$val = [split("", $val)];
		$name = "$gene.$name";
		print "Line error skipped: $line\n" and next if not defined($name) or length($val) == 0;
	#	print "\n\n>$name\n";	
		my $count = 0; my $beg = 0; my $end = 0;
		my $check;
		$check = 1 if $name eq "CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/101566/ccs";
		#next if not defined $check;
		my ($peak, $peakCH, $peakCG, $peakGH, $peakGC) = getPeak($val, $seq{$gene}, $window, $threshold, $check);
		my $t1 = Benchmark->new();
		my ($td) = timestr(timediff($t1, $t0)) =~ /(\-?\d+\.?\d*) wallclock/;
		if ($td > 5) {
			my $persec = int(($linecount-$l0) / $td*10+0.5)/10;
			print STDERR "$peakFilename: Done $LGN$linecount$N in $LCY$persec$N lines per second; ";
			print $outLog "$peakFilename: Done $LGN$linecount$N in $LCY$persec$N lines per second; ";
			$t0 = Benchmark->new();
			$l0 = $linecount;
			print $outLog "$LGN$linecount$N:";
			print STDERR "$LGN$linecount$N:"; my $to = 0;
			foreach my $key (sort keys %pk) {
				next if $key =~ /NO$/;
				my $peak = defined $pk{$key} ? $pk{$key} : 0;
				$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
				print $outLog "$key=$LGN$peak$N, ";
				print STDERR "$key=$LGN$peak$N, ";
			}
			print $outLog "total=$LCY$to$N\n";
			print STDERR "total=$LCY$to$N\n";
		}
		my ($CH, $CG, $GH, $GC) = ("","","","");
		for (my $i = 0; $i < @{$seq{$gene}{seq}}; $i++) {
			$CH .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{CH}[$i] and $peak->{CH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CH}[$i] : "\t1";
			$CG .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{CG}[$i] and $peak->{CG}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CG}[$i] : "\t1";
			$GH .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{GH}[$i] and $peak->{GH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GH}[$i] : "\t1";
			$GC .= $val->[$i] =~ /^[\-\_]$/ ? "\t0" : (defined $peak->{GC}[$i] and $peak->{GC}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GC}[$i] : "\t1";
#			$CH .= $val->[$i] =~ /^[\-\_]$/ ? "\t6" : (defined $peak->{CH}[$i] and $peak->{CH}[$i] =~ /^[\-\_]$/) ? "\t6" : (defined $peak->{CH}[$i] and $peak->{CH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CH}[$i] : "\t7";
#			$CG .= $val->[$i] =~ /^[\-\_]$/ ? "\t6" : (defined $peak->{CG}[$i] and $peak->{CG}[$i] =~ /^[\-\_]$/) ? "\t6" : (defined $peak->{CG}[$i] and $peak->{CG}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{CG}[$i] : "\t7";
#			$GH .= $val->[$i] =~ /^[\-\_]$/ ? "\t6" : (defined $peak->{GH}[$i] and $peak->{GH}[$i] =~ /^[\-\_]$/) ? "\t6" : (defined $peak->{GH}[$i] and $peak->{GH}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GH}[$i] : "\t7";
#			$GC .= $val->[$i] =~ /^[\-\_]$/ ? "\t6" : (defined $peak->{GC}[$i] and $peak->{GC}[$i] =~ /^[\-\_]$/) ? "\t6" : (defined $peak->{GC}[$i] and $peak->{GC}[$i] !~ /^[\! ]$/) ? "\t" . $peak->{GC}[$i] : "\t7";
# 2 is now 2
# 6 is now 6
# bad region is now 6
		}		
		$CH = $peakCH != 0 ? "$name\tPEAK\t$gene\tCH\t$strand\t$CH\n" : "$name\tNOPK\t$gene\tCH\t$strand\t$CH\n";
		$CG = $peakCG != 0 ? "$name\tPEAK\t$gene\tCG\t$strand\t$CG\n" : "$name\tNOPK\t$gene\tCG\t$strand\t$CG\n";
		$GH = $peakGH != 0 ? "$name\tPEAK\t$gene\tGH\t$strand\t$GH\n" : "$name\tNOPK\t$gene\tGH\t$strand\t$GH\n";
		$GC = $peakGC != 0 ? "$name\tPEAK\t$gene\tGC\t$strand\t$GC\n" : "$name\tNOPK\t$gene\tGC\t$strand\t$GC\n";
		$peakCH != 0 ? ($pk{CH} ++ and print $outPEAKCH $CH) : ($pk{CHNO} ++ and print $outNOPKCH $CH);
		$peakCG != 0 ? ($pk{CG} ++ and print $outPEAKCG $CG) : ($pk{CGNO} ++ and print $outNOPKCG $CG);
		$peakGH != 0 ? ($pk{GH} ++ and print $outPEAKGH $GH) : ($pk{GHNO} ++ and print $outNOPKGH $GH);
		$peakGC != 0 ? ($pk{GC} ++ and print $outPEAKGC $GC) : ($pk{GCNO} ++ and print $outNOPKGC $GC);
		#last if $linecount > 100;##	exit 0 if defined $check and $check == 1;
	}
	my $t1 = Benchmark->new();
	my ($td) = timestr(timediff($t1, $t0)) =~ /(\-?\d+\.?\d*) wallclock/;
	my $persec = int(($linecount-$l0) / $td*10+0.5)/10;
	print STDERR "$peakFilename: Done $LGN$linecount$N in $LCY$persec$N lines per second; ";
	print $outLog "$peakFilename: Done $LGN$linecount$N in $LCY$persec$N lines per second; ";
	my $to = 0;
	foreach my $key (sort keys %pk) {
		next if $key =~ /NO$/;
		my $peak = defined $pk{$key} ? $pk{$key} : 0;
		$to = defined $pk{$key . "NO"} ? $pk{$key . "NO"} + $peak : $peak;
		print $outLog "$key=$LGN$peak$N, ";
		print STDERR "$key=$LGN$peak$N, ";
	}
	print $outLog "total=$LCY$to$N\n";
	print STDERR "total=$LCY$to$N\n";
#	print "PEAK: CH=$LGN$peakCH$N, CG=$LGN$peakCG$N, GH=$LCY$peakGH$N, GC=$LCY$peakGC$N\n";
#	print "NOPK: CH=$LGN$peakCH$N, CG=$LGN$peakCG$N, GH=$LCY$peakGH$N, GC=$LCY$peakGC$N\n";
#	exit 0;
	my $peakFilez = "$peakFolder/$peakFilename\_$window\_$threshold\_CG.PEAK";
	print $outLog "footLoop_addition.pl $peakFilez $faFile $gene\n";
	print STDERR "footLoop_addition.pl $peakFilez $faFile $gene\n";
	footPeakAddon::main(($peakFilez, $faFile, $gene, $minDis, $opt_l, $SEQ));
#	system("footLoop_addition.pl $peakFilez $faFile $gene") == 0 or print $outLog "Failed to footLoop_addition.pl $peakFilez $faFile: $!\n";
	next;

=comment
		for (my $i = 0; $i < @val; $i++) {
		#	next if $val[$i] eq "";
		#	print "DIED\n\n$line\n\nVal:" . join(" ", @val[@val-10..@val-1]) . "\n" and die if $val[$i] eq "";
			# 1. Not currently at peak, and current is peak (9)
			if ($peak == 0 and $val[$i] == 9) {
				$peak = 1;
				$totalPeak ++;#= (keys %data);
				my $base = $seq{$gene}[$i]; print $outLog "Died i = $i name=$name total=${@val}\n" and die if not defined($base);
				my $index = int($i / $groupsize);
				$data{$index}{$totalPeak}{beg} = $i;
				$data{$index}{$totalPeak}{name} = $name;
				$bad{$index}{$totalPeak}{beg} = 1 if $i <= 50;
				my $j0 = $i < 500 ? 0 : $i - 500;
				for (my $j = $j0; $j < $i; $j++) {
					last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
					$bad{$index}{$totalPeak}{beg} = 1 if $val[$j] == 6;
					last if $val[$j] == 6;
				}
				push(@{$data{$index}{$totalPeak}{seq}}, $base);
				$count = 0;
	##			print "Peak $totalPeak: name $name: $base";
				$beg = $i; $end = $i;
				$peakcount ++;
			}
			# 2. Currently is peak (peak == 1) but value isn't 9 and total non-peak is equal to maximum buffer length, then peak is end
			elsif ($peak == 1 and (($val[$i] != 9 and $count == $min) or $i == $total - 1 or $i == @val - 1)) {
				my $base = $seq{$gene}[$i]; print $outLog "Died i = $i name=$name total=${@val}\n" and die if not defined($base);
				my $index = int($beg / $groupsize);
				#$end = $i - $min;
				my $indexend = int($end/ $groupsize);
				push(@{$data{$index}{$totalPeak}{seq}}, $base);
				$bad{$index}{$totalPeak}{end} = 1 if $val[$end] == 6 or $i == $total - 1 or $i == @val - 1;
				$data{$index}{$totalPeak}{end} = $end;
	#			for (my $j = $end; $j < $j + 500; $j++) {
				for (my $j = $end; $j < $end + 500; $j++) {
					last if $j == @val - 1 or $j == $total - 1;
					last if $val[$j] == 0 or $val[$j] == 1 or $val[$j] == 2 or $val[$j] == 9;
					$bad{$index}{$totalPeak}{end} = 1 if $val[$j] == 6;
					last if $val[$j] == 6;
				}
				$end{$totalPeak} = $indexend;
	##			print "\ni=$i indexBeg=$index indexEnd=$indexend Peak=$totalPeak name $name: BEG=$beg, END=$end\n\n";
				$count = 0;
				$peak = 0;
				$beg = 0; $end = 0;
			}
			# 3. Currently is peak (peak == 1), but total non-peak length is less than required maximum non-peak length ($min)
			elsif ($peak == 1 and $count < $min) {
				my $base = $seq{$gene}[$i]; 
				print $outLog "\nUndefined base at i=$i gene=$gene, total = $total\n\n" and die if not defined $base;
				$count ++ if $val[$i] != 9;
				$count = 0 if $val[$i] == 9;
				$end = $i if $val[$i] == 9;
				my $index = int($beg / $groupsize);
				push(@{$data{$index}{$totalPeak}{seq}}, $base);
	##			print "$base";
			}
		#	print "\nBEG=$beg, END=$end\n" if $i == @val - 1 or $i == $total - 1;
			last if $i == $total - 1;
		}
		$linecount ++;
		last if $linecount > 10;
#	   last if $totalPeak > 2;
	}
	close $in1;
	exit 0;
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
				my $beg0 = $beg - $minDis < 0 ? 0 : $beg - $minDis;
				my $beg1 = $beg + $minDis > @seq-1 ? @seq-1 : $beg+$minDis;
				my $end0 = $end - $minDis < 0 ? 0 : $end - $minDis;
				my $end1 = $end + $minDis > @seq-1 ? @seq-1 : $end+$minDis;
				my $begSeq = join("", @seq[$beg0..$beg1]);
				my $endSeq = join("", @seq[$end0..$end1]);
				print $outLog "beg0=$beg0; beg1=$beg1; end0=$end0; end1=$end1\n" and die if not defined $seq[$beg0] or not defined($seq[$beg1]);
	#			print "GENE=$gene\tGROUP=$indexBeg.$indexEnd\tBEG=$beg\tEND=$end\tPEAK=$peak\tNAME=$name\tbeg=$begSeq\tend=$endSeq\n";
			}
		}
	}
=cut
	exit 0;
open (my $outLen, ">", "$outDir/$peakFilename.len") or die "Cannot write to $outDir/$peakFilename.len: $!\n";
open (my $outOriginal, ">", "$outDir/$peakFilename.original") or die "Cannot write to $outDir/$peakFilename.original: $!\n";
#my $strand = $peakFile =~ /Pos/i ? "pos" : $peakFile =~ /Neg/i ? "neg" : die "Cannot determine strand for file $peakFile\n";
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
		my $peakEnd = $origBeg + $end + 1; print $outLog "$peak: start ($peakBeg) is less than end ($peakEnd)\n" and die if $peakEnd < $peakBeg;
		#die "Gene = $gene, name = $name, CHR=$origChr, BEG=$origBeg, END=$origEnd, beg = $beg, end = $end, peakBeg = $peakBeg, peakEnd = $peakEnd\n";
		my $len = $end - $beg + 1;
		my $newstrand = $strand eq "pos" ? "+" : $strand eq "neg" ? "-" : ".";
		print $outLen "$gene\t$beg\t$end\t$name\t$len\t$newstrand\n";
		print $outOriginal "$origChr\t$peakBeg\t$peakEnd\t$name\t$len\t$newstrand\n";
		my $badbeg  = $final{$indexBeg}{$indexEnd}{$peak}{badbeg};
		my $badend  = $final{$indexBeg}{$indexEnd}{$peak}{badend};
		my $beg0 = $beg - $minDis < 0 ? 0 : $beg - $minDis;
		my $beg1 = $beg + $minDis > @seq-1 ? @seq-1 : $beg+$minDis;
		my $end0 = $end - $minDis < 0 ? 0 : $end - $minDis;
		my $end1 = $end + $minDis > @seq-1 ? @seq-1 : $end+$minDis;
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
		%kmer = %{kmer($begSeq, "beg", \%kmer)};
		%kmer = %{kmer($endSeq, "end", \%kmer)};
		%kmer = %{kmer($midSeq, "mimd", \%kmer)};
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
print "\t- Relative coordinate bed file: $LGN$outDir/$peakFilename.len$N\n";
print "\t- UCSC Original coordinate bed file: $LGN$outDir/$peakFilename.original$N\n";
print "\t- Kmer odds ratio file: $LGN$outDir/$peakFilename.kmer$N\n";
#print "\n\n--------------------\nKMER:";
#print "\n\n--------------------\n";

open (my $outKmer, ">", "$outDir/$peakFilename.kmer") or die "Cannot write to $outDir/$peakFilename.kmer: $!\n";
my @seqs = @{make_kmer($Kmerz)};
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
		$obs .= "\t$perc";
		$exp .= "\t$percE" if $type eq "beg";
		$odd .= "\t$ratio";
	}
	$obs .= "\n";
	$exp .= "\n" if $type eq "beg";
	$odd .= "\n";
}

print $outKmer ">Observed:\n$obs\n\n";
print $outKmer ">Expected:\n$exp\n\n";
print $outKmer ">OddRatio:\n$odd\n\n";
}
exit 0;

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
	my ($nuc, $seq, $window, $threshold, $check) = @_;
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
	@{$peak->{CH}} = ("!") x (@{$seq->{seq}});
	@{$peak->{CG}} = ("!") x (@{$seq->{seq}});
	@{$peak->{GH}} = ("!") x (@{$seq->{seq}});
	@{$peak->{GC}} = ("!") x (@{$seq->{seq}});
	my $totalseqpos = (keys %{$seq->{loc}{pos1}});
	my $totalseqneg = (keys %{$seq->{loc}{neg1}});
	my $max = $totalseqpos > $totalseqneg ? $totalseqpos : $totalseqneg;
	# initialize
#	return($peak,0,0,0,0); #debug
	for (my $i = 0; $i < $max; $i++) {
		my ($jmin, $jmax) = $i == 0 ? (0, $window) : ($i, $i+1);
		for (my $j = $jmin; $j < $jmax; $j++) {
			if ($j < $totalseqpos - $window) {
				if ($i == 0) {
					my $nucpos1 = $nuc->[$seq->{loc}{pos1}{$j}];
					$con = det_CG_con($con, $nucpos1, 1);
					print $outLog "i=$i, j=$j: nuc.i=$LCY$nucpos1$N, nuc.j(nucpos1)=$nucpos1: file.orig value isn't one of the listed at $LGN" . __LINE__ . "$N\n" if $nucpos1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
					$peak->{CH}[$seq->{loc}{pos1}{$j}] = det_peak($peak->{CH}[$seq->{loc}{pos1}{$j}], $nucpos1, "CH");
					$peak->{CG}[$seq->{loc}{pos1}{$j}] = det_peak($peak->{CG}[$seq->{loc}{pos1}{$j}], $nucpos1, "CG");
				}
				else {
					my ($nucpos0, $nucpos2) = ($nuc->[$seq->{loc}{pos1}{$j-1}], $nuc->[$seq->{loc}{pos1}{$j+$window-1}]);
					$con = det_CG_con($con, $nucpos0, -1, $nucpos2, 1);
					$peak->{CH}[$seq->{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CH}[$seq->{loc}{pos1}{$j+$window-1}], $nucpos2, "CH");
					$peak->{CG}[$seq->{loc}{pos1}{$j+$window-1}] = det_peak($peak->{CG}[$seq->{loc}{pos1}{$j+$window-1}], $nucpos2, "CG");
				}
#				$chunkC .= $nucpos1;
			}
			if ($j < $totalseqneg - $window) {
				if ($i == 0) {
					my $nucneg1 = $nuc->[$seq->{loc}{neg1}{$j}];
					$con = det_CG_con($con, $nucneg1, 1);
					print $outLog "i=$i, j=$j: nuc.i=$LCY$nucneg1$N, nuc.j(nucneg1)=$nucneg1: file.orig value isn't one of the listed at $LGN" . __LINE__ . "$N\n" if $nucneg1 !~ /^[BCDcdMNFEeOJGHghUVKIiW\-_\.TA123PQR456XYZ]$/;
					$peak->{GH}[$seq->{loc}{neg1}{$j}] = det_peak($peak->{GH}[$seq->{loc}{neg1}{$j}], $nucneg1, "GH");
					$peak->{GC}[$seq->{loc}{neg1}{$j}] = det_peak($peak->{GC}[$seq->{loc}{neg1}{$j}], $nucneg1, "GC");
				}
				else {
					my ($nucneg0, $nucneg2) = ($nuc->[$seq->{loc}{neg1}{$j-1}], $nuc->[$seq->{loc}{neg1}{$j+$window-1}]);
					$con = det_CG_con($con, $nucneg0, -1, $nucneg2, 1);
					$peak->{GH}[$seq->{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GH}[$seq->{loc}{neg1}{$j+$window-1}], $nucneg2, "GH");
					$peak->{GC}[$seq->{loc}{neg1}{$j+$window-1}] = det_peak($peak->{GC}[$seq->{loc}{neg1}{$j+$window-1}], $nucneg2, "GC");
				}
#				$chunkG .= $nucneg1;
			}
#=cut
		}
		if ($i == 0 and defined $check) {
#			print "i=$i, " . join("", @{$peak->{CH}}[0..$window-1]) . "\n";
		}
		my ($CG, $CH, $GC, $GH) = (0,0,0,0);
		($peak, $CG, $CH, $GC, $GH) = isPeak($i, $con, $peak, $window, $threshold, $seq, $totalseqpos, $totalseqneg, $check);
		$peakCH += $CH;
		$peakCG += $CG;
		$peakGC += $GC;
		$peakGH += $GH;
	}
	return ($peak, $peakCH, $peakCH, $peakGC, $peakGH);
}
sub det_peak { # time= 10 line/s per strand (49 l/s -> 39 l/s)
	my ($peak, $nuc, $type) = @_;
#	return 1; #debug
	return $peak if $peak ne '!';
	my @N =	$type eq 'CH' ? ('0', '0'   , 'cd', 'CDMNPQ', '12' , 'B-' ) :
				$type eq 'CG' ? ('e', 'EQR' , 'cd', 'CDMNPQ', '123', 'BF-') :
				$type eq 'GH' ? ('0', '0'   , 'gh', 'GHUVXY', '45' , 'J-' ) :
				$type eq 'GC' ? ('i', 'IWZ' , 'gh', 'GHUVXY', '456', 'JK-') : 
				die "det_peak TYPE IS NOT CH/CG/GH/GC? ($type) at line $LGN" . __LINE__ . "$N\n";
	return 7 if $nuc eq $N[0]; # conv CG
	return 5 if $nuc eq $N[1]; # nonconverted CG
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
sub kmer {
	my ($seq, $type, $kmer) = @_;
	my %kmer = %{$kmer};
	$seq = uc($seq);
	my @seqz = split("", $seq);
	for (my $i = 0; $i < @seqz-($Kmerz-1); $i++) {
		my $kmer;
		for (my $j = $i; $j < $i+$Kmerz; $j++) {
			$kmer .= $seqz[$j];
		}
		$kmer{$type}{$kmer} ++;
		$kmer{$type}{total} ++;
	}
#	for (my $i = 0; $i < @seqz-1; $i++) {
#		my $kmer = "$seqz[$i]$seqz[$i+1]";
#		print "$type $kmer = $kmer{$type}{$kmer}\n";
#	}
	return(\%kmer);
}
sub record_options {
   my ($opts, $logs, $other, $outLog, $logFile, $date, $uuid) = @_;
   my $optPrint = "$0";
	my $optShort = "$0";
   foreach my $opt (sort keys %{$opts}) {
		next if $opt eq "o";
      my $val = $opts->{$opt};
      $optPrint .= $val eq "MYTRUE" ? " -$opt" : " -$opt $val";
		$optShort .= defined ($logs->{$opt}) ? "" : $val eq "MYTRUE" ? " -$opt" : " -$opt $val";
   }
   my  $param = "

$YW<-----------------------------------------------$N
${YW}Initializing...$N
   
>Run Params:
Date                $LCY: $date; $N
Run ID              $LCY: $uuid; $N
Run script full     $LCY: $optPrint; $N
Run script short    $LCY: $optShort; $N
>Options:
-d minDis  $LCY: $opts->{d}; $N
-s grpSize $LCY: $opts->{s}; $N
-k Dist    $LCY: $opts->{k}; $N
-K Kmer    $LCY: $opts->{K}; $N
-t thrshld $LCY: $opts->{t}; $N
-w window  $LCY: $opts->{w}; $N

>Run Params from footLoop.pl logfile=$logFile
footLoop Date       $LGN: $other->{date}; $N
footLoop Run ID     $LGN: $other->{uuid}; $N
footLoop Run script $LGN: $other->{runscript}; $N
footLoop md5        $LGN: $other->{md5}; $N
footLoop origDir    $LGN: $other->{origDir}; $N
>Options from footLoop.pl logfile=$logFile\n";
	 my $optcount = 0;
	foreach my $opt (sort keys %{$opts}) {
		next if not defined $logs->{$opt}; next if $opt eq "origDir";
		$optcount ++;
		$param .= "$optcount. -$opt $LGN: $opts->{$opt}; $N\n";
	}
	$param .= "\n\n$YW----------------------------------------------->$N\n\n";

   print $outLog $param;
   print STDERR $param;
}

sub parse_footLoop_logFile {
	my ($logFile, $date, $uuid) = @_;
	my @line = `cat $logFile`;

	my $defOpts = set_default_opts();
	
	# %others contain other parameters from footLoop that aren't options (e.g. uuid)
	# %log TBA

	my ($other, $log); 
	foreach my $line (@line[0..@line-1]) {
		if ($line =~ /^Date\s*:/) {
			($other->{date}) = $line =~ /^Date\s+:\s*([a-zA-Z0-9\:\-].+)$/;
			print "Undefined date from $line\n" and die if not defined $other->{date};
		}
		if ($line =~ /^Run ID\s*:/) {
			($other->{uuid}) = $line =~ /^Run ID[ \t]+:[ \t]+([a-zA-Z0-9]+.+)$/;
			print "Undefined uuid from $line\n" and die if not defined $other->{uuid};
		}
		if ($line =~ /^Run script\s*:/) {
			($other->{runscript}) = $line =~ /^Run script[ ]+:(.+)$/;
			print "Undefined runscript from $line\n" and die if not defined $other->{runscript};
			$defOpts = parse_runscript($defOpts, $other->{runscript});
		}
	}
	my ($outDirs) = `tail -n 10 $logFile | grep 'Output: '`; chomp($outDirs);
	my @outDirs = split("\/", $outDirs);
	my $outDir = ""; my $temp = "";
	for (my $i = 0; $i < @outDirs; $i++) {
		my $temp = $outDir . "/" . $outDirs[$i];
		$outDir = -d $temp ? $temp : $outDir;
	}
	$other->{origDir} = $outDir;
	$defOpts->{origDir} = $outDir;
	($other->{md5}) = $outDir =~ /\.0_orig_([a-z0-9]+)/i; print $outLog "Can't find outdir md5 (outdir=$defOpts->{o})\n" and die if not defined $other->{md5};
#	die "opt i = $defOpts->{i} others=$other->{runscript}\n";
	my ($indexFilename) = getFilename($defOpts->{i}, "full");
	my ($faFiles) = `grep "Sequence has been parsed from fasta file" $logFile` =~ /^.+fasta file (.+)$/;

	my @faFile = split("\/", $faFiles);
	my $faFile = ""; $temp = "";
	for (my $i = 0; $i < @faFile; $i++) {
		$temp = $faFile . "/" . $faFile[$i];
		$faFile = -d $temp ? $temp : $faFile;
	}
	$faFile .= "/$indexFilename\_$defOpts->{x}\_$defOpts->{y}\_bp.bed.fa";
	$defOpts->{g} = $faFile;
	print $outLog "Fasta file $faFile does not exist!\n" and die if not -e $faFile;
	open (my $outLog, ">", "$outDir/footPeak_logFile.txt") or print "Failed to write to $outDir/footPeak_logFile.txt: $!\n" and exit 1;
	record_options($defOpts, $log, $other, $outLog, $logFile);
	return($defOpts, $outLog);
}

sub parse_runscript {
	my ($defOpts, $runscripts) = @_;
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
				$val = "TRUE";
			}
			if (not defined $defOpts->{$opt}) {
				print "parse_runscript: opt $opt does not exist!\n" and next;
			}
		}
		elsif ($i == 0) {
			print "script = $runscripts[0]\n" and next;
		}
		else {
			print "parse_runscript: $runscripts[$i] nexted\n" and next;
		}
		if ($val =~ /^\-[a-zA-Z]$/) {
			my ($val1) = $val =~ /^\-(.)$/;
			if (defined $defOpts->{$val1}) {
				$val = "TRUE";
			}
			else {
				$i ++;
			}
		}
		else {
			$i ++;
		}
		$defOpts->{$opt} = $val;
		$log->{$opt} = 1;
	}
	my $defOptsCount = 0;
	my $defOptsPrint = "\%defOpts = (";
	foreach my $opt (sort keys %{$defOpts}) {
		#$defOptsPrint .= "\'$opt\' => \'$defOpts->{$opt}\',";
		$defOptsPrint .= "\n\t\t" if $defOptsCount % 4 == 0;
		$defOptsPrint .= "\'$opt\' => \'\$opt\_$opt\',\t" if defined $log->{$opt};
		$defOptsPrint .= "\'$opt\' => \'$defOpts->{$opt}\',\t" if not defined $log->{$opt};
		$defOptsCount ++;
	}
	$defOptsPrint =~ s/,\t$/\n/;
	$defOptsPrint .= ");";
	print "$defOptsPrint\n";
	die;
	for (my $i = 1; $i < @runscripts; $i+=2) {
		my $defOpts = $runscripts[$i];
		my ($opt, $val) = ($runscripts[$i], $runscripts[$i+1]);
		if ($opt =~ /^\-[a-zA-Z0-9]$/ and $val =~ /^\-[a-zA-Z0-9]$/ and defined $runscripts[$i+2] and $runscripts[$i+2] !~ /^\-[0-9]$/) {
			($opt, $val) = ($opt, "TRUE");
			$i -= 1;
		}
		($opt) = $defOpts =~ /^(.)$/ if not defined $opt;
		$opt =~ s/^\-//;
		next if not defined $opt;
		$defOpts->{$opt} = $val;
		#$log->{$opt} = $val;
		print "-$opt $val\n";
	}
	return $defOpts;
}

sub set_default_opts {
	# To define default values as well as record purposes, we create 4 hashes, %defOpts, %usrOpts, %other, and %log.
	# %defOpts contains options that came from footLoop logfile.txt and came from user inputs from this script (footPeak.pl)
	# %defOpts stores default values, %usrOpts stores user inputs.
	my %defOpts = (
					'x' => 0,	'y' => 0,	'g' => '',	'i' => '',	'q' => 0, 	'r' => '', 
					'L' => 0,	't' => 65,	'w' => 20,	'd' => 250, 's' => 200, 'k' => 50,
					'K' => 2, 	'n' => '', 	'l' => ''
					);
		
	# %usrOpts contains options only from this script (footPeak.pl) 
	# For each undefined values in %usrOpts (user didn't input anything), we use default value in %defOpts.
	my %usrOpts = (
					't' => $opt_t, 'w' => $opt_w, 'd' => $opt_d, 's' => $opt_s, 
					'k' => $opt_k, 'K' => $opt_K, 'n' => $opt_n, 'l' => $opt_l
					);
	# set opt in footPeak's %defOpts to be what user inputted based on %usrOpts.
	# if doesn't exist, then use default value in %defOpts.
	foreach my $opt (sort keys %usrOpts) {
		$defOpts{$opt} = $usrOpts{$opt} if defined $usrOpts{$opt};
	}
	return \%defOpts;
}

sub parse_index {
	my ($indexFile, $outDir, $leftBuf, $riteBuf) = @_;
	my ($folder, $fileName) = getFilename($indexFile, "folder");
	my $indexFile2 = $outDir . "/$fileName\_$leftBuf\_$riteBuf\_bp.bed";
	return $indexFile2;
	#system("bedtools_bed_change.pl -x $leftBuf -y $riteBuf -i $indexFile -o $indexFile\_$leftBuf\_$riteBuf.bed > /dev/null 2>&1");

#	my @LINE = `cat $indexFile\_$leftBuf\_$riteBuf.bed`;
#	print "\n$YW Parsing index file $indexFile\_$leftBuf\_$riteBuf.bed$N\n";
#	foreach my $line (@LINE) {
#		chomp($line);
#		my ($chr, $beg, $end, $name) = split("\t", $line);
#		print "\tParsed:Gene=$CY$name$N, COORD=$CY$chr\t$beg\t$end\t$name$N\n\n";
#		$coor{uc($name)}{chr} = $chr;
#		$coor{uc($name)}{beg} = $beg;
#		$coor{uc($name)}{end} = $end;
#	}
#	system("fastaFromBed -fi $faFile -bed $indexFile\_$leftBuf\_$riteBuf.bed -fo $indexFile\_$leftBuf\_$riteBuf.fa -name");
#	$faFile = "$indexFile\_$leftBuf\_$riteBuf.fa";
#	return($faFile);
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
	#print "\nK=$Kmerz; Kmer = " . join(",", @preword[0..$max]);
	#print "\n" if @preword <= 9;
	#print "..." . "(total =$LGN " . scalar(@preword) . "$N)\n\n" if @preword > 9;
	return(\@preword);
}

sub check_sanity {
	my ($footFolder) = @_;


my $usage = "
Usage: $YW$0$N -i <footLoop output folder>
#-g$CY <genomic fasta>$N -i$LPR <UNMODIFIED geneIndexes.bed>$N -p$LGN <Peak file>$N -n$YW <Output Folder>$N

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
-d: Maxmimum distance between 2 peaks in basepair [100]
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

";

print $usage . $usage_long and die if defined $opt_h;
print $usage and die if not defined $footFolder or not -d $footFolder;
print $outLog "Please run footloop.pl first! ($footFolder/logFile.txt does not exists)\n" and die if not -e "$footFolder/logFile.txt";

return $usage;
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
__END__
my ($faFile, $indexFile, $peakFile, $x, $y, $min, $groupsize, $minDis, $outDir) = ($opt_g, $opt_i, $opt_p, $opt_x, $opt_y, $opt_d, $opt_s, $opt_k, $opt_n)

__END__
		print "\n$name\tCG\t$strand\t";		
		foreach my $loc (sort {$a <=> $b} @{$data->{CG}}) {
			print "$data->{CG}{$loc}";
		}		
		print "\n$name\tGH\t$strand\t";		
		foreach my $loc (sort {$a <=> $b} @{$data->{GH}}) {
			print "$data->{GH}{$loc}";
		}		
		print "\n$name\tGC\t$strand\t";		
		foreach my $loc (sort {$a <=> $b} @{$data->{GC}}) {
			print "$data->{GC}{$loc}";
		}		
		print "\n";

