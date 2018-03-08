#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';
use vars qw($opt_v $opt_n);
getopts("n:");

my ($footPeakFolder) = ($opt_n);
die "\nusage: $YW$0$N -n $CY<footPeakFolder>$N\n\n" unless defined $opt_n and -d $opt_n;

my ($faFile) = <$footPeakFolder/*.fa>;
die if not defined $faFile or not -e $faFile;
my %fa;
open (my $faIn, "<", $faFile) or die;
my $fasta = new FAlite($faIn);
while (my $entry = $fasta->nextEntry()) {
	my $def = $entry->def(); $def =~ s/^>//;
	my $seq = $entry->seq();
	$fa{$def} = $seq;
}

my @bed = <$footPeakFolder/FOOTCLUST/.TEMP/*.bed>;
makedir("$footPeakFolder/FOOTCLUST/.DREME/");
open (my $outz, ">", "$footPeakFolder/FOOTCLUST/dreme.txt");
foreach my $bedFile (sort @bed) {
	my ($folder1, $fileName1) = getFilename($bedFile, "folderfull");
	my ($gene, $strand, $window, $thres, $type) = $fileName1 =~ /^(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC)\./;
	print "File=$fileName1, gene=$gene, strand=$strand, window=$window, thres=$thres, type=$type\n";
	open (my $in, "<", $bedFile) or die;
	my $linecount = 0;
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /^#/;
		print "Done $linecount\n" if $linecount % 5 == 0;
		$linecount ++;
		my($chr, $beg, $end, $cluster, $total_peak, $strand2) = split("\t", $line);
		my $seq = $fa{$chr}; die if not defined $fa{$chr};
		my $len = length($seq);
		if ($end - $beg < 20) {
			my $mid = int(($beg + $end) / 2+0.5);
			$beg = $mid - 10;
			$end = $mid + 10;
		}
		if ($end - $beg > 50 and $cluster !~ /WHOLE/) {
			my $mid = int(($beg + $end) / 2+0.5);
			$beg = $mid - 25;
			$end = $mid + 25;
		}
		my $origBed = "$footPeakFolder/FOOTCLUST/.DREME/$gene\_$strand\_$window\_$thres\_$type\_$cluster.orig";
		my $shufBed = "$footPeakFolder/FOOTCLUST/.DREME/$gene\_$strand\_$window\_$thres\_$type\_$cluster.shuf";
		my $dremeOut = "$footPeakFolder/FOOTCLUST/.DREME/$gene\_$strand\_$window\_$thres\_$type\_$cluster/";
		open (my $outorig, ">", $origBed) or die;
		print $outorig "$chr\t$beg\t$end\t$cluster.orig\t$total_peak\t$strand\n";
		my ($num, $clust) = $cluster =~ /^(\d+).(\w+)$/;
		close $outorig;
		open (my $outshuf, ">", $shufBed) or die;
		my $lenz = $len - $end;
		my ($seq1) = $seq =~ /^(.{$beg})/;
		my ($seq2) = $seq =~ /(.{$lenz})$/;
		my $len0 = $end - $beg;
		die "undefined seq2 len=$len beg=$beg end=$end\n" if not defined $seq2;
		my ($len1, $len2) = (length($seq1), length($seq2));
		for (my $i = 0; $i < 100; $i++) {
			my $rand = rand($len);
			my $currseq = ($len1 < $len0 and $len2 < $len0 and $len2 < $len1) ? $seq1 :
							  ($len1 < $len0 and $len2 < $len0 and $len2 > $len1) ? $seq2 :
								($len1 > $len0 and $len2 < $len2) ? $seq1 :
								($len1 < $len0 and $len2 > $len2) ? $seq2 :
								$rand > $len1 ? $seq2 : $seq1;
			my $currbeg = ($len1 < $len0 and $len2 < $len0) ? 1 : int(rand(length($currseq)-$len0)); 
			$currbeg = 1 if $currbeg < 0;
			my $currend = ($len1 < $len0 and $len2 < $len0) ? length($currseq) : $currbeg + $len0;
			print $outshuf "$chr\t$currbeg\t$currend\t$cluster.shuf.1\t$total_peak\t$strand\n";
		}
#		print $outshuf "$chr\t0\t$beg\t$cluster.shuf.1\t$total_peak\t$strand\n";
#		print $outshuf "$chr\t$end\t$len\t$cluster.shuf.2\t$total_peak\t$strand\n";
		system("fastaFromBed -fi $faFile -bed $origBed -fo $origBed.fa -s -name");
		system("fastaFromBed -fi $faFile -bed $shufBed -fo $shufBed.fa -s -name");
		my $dreme = `dreme -p $origBed.fa -n $shufBed.fa -oc $dremeOut | grep \"Best RE\"`;
		chomp($dreme);
		my ($motif, $pval, $eval) = $dreme =~ /^Best RE was (.+) \w+ p-value= (.+) E-value= (.+) Unerased/;
		print $outz "$chr\t$beg\t$end\t$cluster\t$total_peak\t$gene\t$strand\t$window\t$thres\t$type\t$num\t$clust\t$motif\t$pval\t$eval\n";
	}
	close $in;
}
__END__


my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
}
close $in1;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

