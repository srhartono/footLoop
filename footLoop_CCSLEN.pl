#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v);
getopts("v");

my ($input1) = @ARGV;
die "\nusage: $YW$0$N $CY<input1>$N\n\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");

my $max = 0;
my %data;
my $in1;
open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.gz$/;
open ($in1, "zcat $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.gz$/;
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
	$line = <$in1>; chomp($line);
	$data{int(length($line)/1000)} ++;
	$max = int(length($line)/1000) if $max < int(length($line)/1000);
	$line = <$in1>;
	$line = <$in1>;
}
close $in1;
open (my $out1, ">", "$input1.ccslen");
foreach my $len (sort {$a <=> $b} keys %data) {
	my $len0 = $len * 1000;
	my $len1 = ($len + 1) * 1000;
	#print "$LGN$len0$N\t$LCY$len1$N\t$LPR$data{$len}$N\n";
	print $out1 "$fileName1\t$len0\t$len1\t$data{$len}\n";
}
close $out1;
