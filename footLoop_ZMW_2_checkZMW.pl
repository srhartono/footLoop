#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_i $opt_b);
getopts("vi:b:");

my ($input1, $bamFile) = ($opt_i, $opt_b);
die "\nusage: $YW$0$N -i $CY<input1.fq.gz>$N -b$LGN <subreads.bam>$N\n\n" unless defined $opt_i and defined $opt_b and -e $opt_i and -e $opt_b;

my %bam;
print "1. doing $LCY$bamFile$N\n";
my $Linecount = 0;
open (my $in0, "samtools view $bamFile|") or die "Cannot read from $bamFile: $!\n";
while (my $line = <$in0>) {
	chomp($line);
	my $print = 1;
	$Linecount ++;
	print "Done $Linecount\n" if $Linecount % 1000000 == 0;
	my @arr = split("\t", $line); next if @arr < 10;
	my ($zmw) = $arr[0] =~ /\/(\d+)\/\d+_\d+/; die "\n\n$line\n\n" if not defined $zmw;
	$bam{$zmw} = 1;
}
my $totalbam = (keys %bam);
print "total bam = $LGN$totalbam$N\n";

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");


my %total;
my %data;
my %line;
open (my $in1, "zcat $input1|") or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$input1.out") or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my $print = 1;
	my ($zmw) = $line =~ /\/(\d+)\/ccs/; die if not defined $zmw;
	if (defined $data{$zmw}) {
		print "$LGN$zmw$N: $line\n";
		$total{zmw}{$zmw} ++;
		$total{total} ++;
		$total{zmwline}{$line} = 1;
		$print = 0;
	}
	$data{$zmw}++;
	$line{$zmw}{$line} ++;
	print "zmw doesnt exist in bam: $zmw\n" if not defined $bam{$zmw};
	$bam{$zmw} = 0;
	my $line2 = <$in1>;
	my $line3 = <$in1>;
	my $line4 = <$in1>;
	if ($print == 1) {
		print $out1 "$line\n$line2$line3$line4";
	}
}
close $in1;
$total{total} = 0 if not defined $total{total};
my $zmwtot = (keys %{$total{zmw}});
my $zmwdata = (keys %data);
foreach my $zmw (sort {$data{$b} <=> $data{$a}||$a cmp $b} keys %data) {
	last if $data{$zmw} == 1;
	foreach my $line (sort keys %{$line{$zmw}}) {
		print "$zmw\t$data{$zmw}\t$line{$zmw}{$line}\t$line\n";
	}
}

print "total = $total{total} (unique = $zmwtot) / $zmwdata\n";
print "total bam = $LGN$totalbam$N\n";

my $notexist = 0;
open (my $out2, ">", "$bamFile.notinccs") or die;
foreach my $zmw (sort {$bam{$b} <=> $bam{$a}} keys %bam) {
	$notexist ++ if $bam{$zmw} == 1;
	print $out2 "$zmw\n" if $bam{$zmw} == 1;
	last if $bam{$zmw} == 0;
}
close $out2;

print "exist in bam, but not in ccs: $notexist\n";
