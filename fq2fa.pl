#!/usr/bin/perl

use strict; use warnings;

my ($input1) = @ARGV;
die "usage: $0 <input1.fq>\n" unless @ARGV;

my $linecount = 0;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$input1.fa") or die "Cannot write to $input1.fa $!\n";
while (my $line = <$in1>) {
	chomp($line);
	$linecount ++;
	print "Done $linecount\n" if $linecount % 1000 == 0;
	# name
	my $name = $line; $name =~ s/^\@//;
	# seq
	$line = <$in1>;chomp($line);
	my $seq  = $line;
	# plus
	$line = <$in1>;chomp($line);
	# plus
	$line = <$in1>;chomp($line);
	my $qual = $line;
	print $out1 ">$name\n$seq\n";
}
close $in1;
close $out1;
