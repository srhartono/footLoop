#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use Thread; use Thread::Queue;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_V);
getopts("vV");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';

my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";
my @version = `cd $footLoopDir && git log | head `;
my $version = "UNKNOWN";
foreach my $line (@version[0..@version-1]) {
   if ($line =~ /^\s+V\d+\.?\d*\w*\s*/) {
      ($version) = $line =~ /^\s+(V\d+\.?\d*\w*)\s*/;
   }
}
if (not defined $version or (defined $version and $version eq "UNKNOWN")) {
   ($version) = `cd $footLoopDir && git log | head -n 1`;
}
if (defined $opt_v) {
   print "\n\n$YW$0 $LGN$version$N\n\n";
   exit;
}
my %undef;
my %coor;
my $barcode = "/home/mitochi/pacbio_indexes/set1and2indexes.bed";
my @line = `cat $barcode`;
foreach my $line (@line) {
	chomp($line);
	next if $line =~ /^#/;
	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	$gene = uc($gene);
	$coor{$gene} = $strand;
	$undef{$gene}{exist} = 1;
}
my ($folder) = @ARGV;
die "\nusage: $YW$0$N [-V verbose; to show all result] $CY<folder containing PEAK and PEAKNEG>$N\n\n" unless @ARGV;

my @peak = -e "$folder/PEAK/" ? <$folder/PEAK/*.bed> : ();
my @peakneg = -e "$folder/PEAKNEG/" ? <$folder/PEAKNEG/*.bed> : ();

foreach my $input1 (sort @peak) {
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	my ($label, $gene, $strand, $window, $thres, $type) = parseName($fileName1);
	my $strand2 = $strand =~ /^(Pos|Neg|Unk)$/ ? $strand : $strand eq "+" ? "Pos" : $strand eq "-" ? "Neg" : "Unk";
	$gene = uc($gene);
	my ($num) = $label =~ /PCB(\d+)$/;
	$undef{$gene}{exist} = 0;
	$undef{$gene}{label}{$num} = 1;
	my $goodstrand = $coor{$gene}; $goodstrand = "FALSE" if not defined $goodstrand;
	my $goodstrand0 = $goodstrand;
	$goodstrand = $goodstrand =~ /^(Pos|Neg|Unk)$/ ? $goodstrand : $goodstrand eq "+" ? "Pos" : $goodstrand eq "-" ? "Neg" : $goodstrand;
	if ($goodstrand ne "FALSE" and $goodstrand ne $strand2) {
		$goodstrand = "${LRD}BAD$N";
		print "$LPR PEAK$N $label\t$gene\t$strand\t$window\t$thres\t$type\t(strand2 $strand2 == $goodstrand? ($goodstrand0))\n";
	} else { $goodstrand = "${LGN}GOOD$N";}
	print "$LPR PEAK$N $label\t$gene\t$strand\t$window\t$thres\t$type\t(strand2 $strand2 == $goodstrand? ($goodstrand0))\n" if defined $opt_V;

	die "Peak $gene Pos not CH\n" if $strand eq "Pos" and $type ne "CH";
	die "Peak $gene Neg not GH\n" if $strand eq "Neg" and $type ne "GH";
	
}

foreach my $input1 (sort @peak) {
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	my ($label, $gene, $strand, $window, $thres, $type) = parseName($fileName1);
	my $strand2 = $strand =~ /^(Pos|Neg|Unk)$/ ? $strand : $strand eq "+" ? "Pos" : $strand eq "-" ? "Neg" : "Unk";
	$gene = uc($gene);
	my $goodstrand = $coor{$gene}; $goodstrand = "FALSE" if not defined $goodstrand;
	my $goodstrand0 = $goodstrand;
	$goodstrand = $goodstrand =~ /^(Pos|Neg|Unk)$/ ? $goodstrand : $goodstrand eq "+" ? "Pos" : $goodstrand eq "-" ? "Neg" : $goodstrand;
	if ($goodstrand ne "FALSE" and ($goodstrand eq "Pos" and $strand2 ne "Neg") and ($goodstrand eq "Neg" and $strand2 ne "Pos")) {
		$goodstrand = "${LRD}BAD$N";
		print "$YW PEAKNEG$N $label\t$gene\t$strand\t$window\t$thres\t$type\t(strand2 $strand2 == $goodstrand? ($goodstrand0))\n";
	} else { $goodstrand = "${LGN}GOOD$N";}
	print "$YW PEAKNEG$N $label\t$gene\t$strand\t$window\t$thres\t$type\t(strand2 $strand2 == $goodstrand? ($goodstrand0))\n" if defined $opt_V;
	die "PeakNeg $gene Pos not CH\n" if $strand eq "Pos" and $type ne "CH";
	die "PeakNeg $gene Neg not GH\n" if $strand eq "Neg" and $type ne "GH";
}

print "gene,total";
for (my $i = 1; $i < 20; $i++) {
	print ",PCB$i";
}
print "\n";
foreach my $gene (sort keys %coor) {
	if ($undef{$gene}{exist} eq 0) {
		my @labels;
		foreach my $num (sort  {$a <=> $b} keys %{$undef{$gene}{label}}) {
			push (@labels, $num);
		}
#		print "$gene," . scalar(@labels) . ",PCB" . join(",PCB", @labels) . "\n";
		print "$gene," . scalar(@labels);
		for (my $i = 1; $i < 20; $i++) {
			defined $undef{$gene}{label}{$i} ? print ",TRUE" : print ",";
		}
		print "\n";
	}
	else {
		print "$gene,0";
		for (my $i = 1; $i < 20; $i++) {
			print ",";
		}
		print "\n";
	}
}

__END__
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");


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

