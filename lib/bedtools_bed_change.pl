#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_x $opt_y $opt_o $opt_a $opt_b $opt_i $opt_s $opt_f $opt_m $opt_c $opt_r $opt_q $opt_v);
getopts("i:x:o:y:abcsf:mrqv");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . "/lib";
   push(@INC, $libPath);
}

use myFootLib;
use FAlite;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);
my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

my $usage = " 

-----------------
$YW $0 $version_small $N
-----------------

Usage: $0 [option] -i bed file

options:
-v: get version
-a: Get start of gene (strand specific) and offset accordingly
-b: Get end of gene (strand specific) and offset accordingly
-c: Get middle of gene and offset accordingly
-s: Disable strand specific (default: on)
-x: start offset from current (default: 0)
-y: end offset from current (default: 0)
-o: output (name)
-f: Filter genes less than ths length (default: take all)
-m: Stop filtering out genes which after subtraction has coordinates less than 0 (Default: Filtered)
-r: randomly add/subtract x/y
-q: Quiet

E.g. you want all bed file
chr1	5000	6000	name	0	-
to become
chr1	4000	10000	name	0	-
then:

$0 -i foo.bed -x -1000 -y 4000 -o bar.bed 

E.g. +/- 1kb region of TSS (strand specific)
$0 -i foo.bed -x -1000 -y 1000 -o bar.bed

";

die $usage unless defined($opt_i);

my $input  = $opt_i;
my $x_off  = defined($opt_x) ? $opt_x : 0;
my $y_off  = defined($opt_y) ? $opt_y : 0;
my $output = defined($opt_o) ? $opt_o : "$input.bed";
my $filter = defined($opt_f) ? $opt_f : -99999999999999;

print "\tOutput = $output\n" if not defined $opt_q;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", $output) or die "Cannot write to $output: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /track/;
	next if $line =~ /\#/;
	my ($chr, $start, $end, $name, $val, $strand, @others) = split("\t", $line);
	next if not defined($chr);
	$val = 0 if not defined($val);
	next if $end - $start + 1 < $filter;

	if (not defined($strand)) {
		$strand = "+";
	}
	else {
		if ($strand ne "-" and $strand ne "+") {
			if (defined($others[0]) and ($others[0] eq "+" or $others[0] eq "-")) {
				my $temp = $strand;
				$strand = $others[0];
				$others[0] = $temp;
			}
			else {
				die "Strand information is incorrect (strand = $strand) at line $line\n";
			}
		}
	}
	my $others = join("\t", @others) if defined($others[0]);
	my ($newstart, $newend);
	if (not defined($opt_a) and not defined($opt_b) and not defined($opt_c)) {
		if (not defined($opt_s)) {
			if ($strand eq "+" or $strand eq "1" or $strand eq "F") {
				$newstart = $start + $x_off;
				$newend   = $end   + $y_off;
			}
			elsif ($strand eq "-" or $strand eq "-1" or $strand eq "R") {
				$newstart = $start - $y_off;
				$newend   = $end   - $x_off;
			}
		}
		else {
			$newstart = $start + $x_off;
			$newend   = $end   + $y_off;
		}
	}
	elsif (defined($opt_a)) {
		my $pos = $strand eq "+" ? $start : $end;
		$newstart = $strand eq "+" ? $pos + $x_off : $pos - $y_off;
		$newend   = $strand eq "+" ? $pos + $y_off : $pos - $x_off;
	}
	elsif (defined($opt_b)) {
		my $pos = $strand eq "+" ? $end : $start;
		$newstart = $strand eq "+" ? $pos + $x_off : $pos - $y_off;
		$newend   = $strand eq "+" ? $pos + $y_off : $pos - $x_off;
	}
	elsif (defined($opt_c)) {
		my $pos = int(($start + $end) / 2);
		$newstart = $strand eq "+" ? $pos + $x_off : $pos - $y_off;
		$newend   = $strand eq "+" ? $pos + $y_off : $pos - $x_off;
	}
	$newstart = 1 if $newstart < 1;
	$newend = 1 if $newend < 1;

	if (($newstart < 1 or $newend < 1) and not $opt_m) {
		print "Skipped [$chr $start $end (new = $newstart $newend)] because start and/or end is less than 1\n";
		next;
	}
	if ($newstart > $newend) {
		#print "Skipped [$chr $start $end] because start is bigger than end!\n";
		my $tempstart = $newstart;
		$newstart = $newend;
		$newend   = $tempstart;
	}
	if (not defined($name)) {
		if ($opt_r) {
			my $posbuffer = (0.1 * ($newend - $newstart) < 10000) ? 10000 : (0.1 * ($newend - $newstart) > 50000) ? 50000 : 0.1 * ($newend - $newstart);
			$posbuffer = rand() <= 0.5 ? -1 * int(rand($posbuffer)) : int(rand($posbuffer));
			my $orignewstart = $newstart; $newstart = $newstart + $posbuffer;
			my $orignewend = $newend; $newend = $newend + $posbuffer;
			if ($newstart <= 0 or $newend <= 0) {
				$newstart = $orignewstart;
				$newend = $orignewend;
			}
		}
		die "UNDEFINED CHR AT LINE $line\n" unless defined($chr);
		print $out "$chr\t$newstart\t$newend\n" if not defined($name) and $newstart < $newend;
		print $out "$chr\t$newend\t$newstart\n" if not defined($name) and $newstart >= $newend;
	}
	else {
		if ($opt_r) {
			my $posbuffer = (0.1 * ($newend - $newstart) < 10000) ? 10000 : (0.1 * ($newend - $newstart) > 50000) ? 50000 : 0.1 * ($newend - $newstart);
			$posbuffer = rand() <= 0.5 ? -1 * int(rand($posbuffer)) : int(rand($posbuffer));
			my $orignewstart = $newstart; $newstart = $newstart + $posbuffer;
			my $orignewend = $newend; $newend = $newend + $posbuffer;
			if ($newstart <= 0 or $newend <= 0) {
				$newstart = $orignewstart;
				$newend = $orignewend;
			}
			#my $posbuffer = (0.1 * ($newend - $newstart) < 10000) ? 10000 : (0.1 * ($newend - $newstart) > 50000) ? 50000 : 0.1 * ($newend - $newstart);
			#my $orignewstart = $newstart;
			#$newstart = rand() <= 0.5 ? (-1 * int(rand($posbuffer))) + $newstart : int(rand($posbuffer)) + $newstart; $newstart = $orignewstart if $newstart <= 0;
			#my $orignewend = $newend;
			#$newend = rand() <= 0.5 ? (-1 * int(rand($posbuffer))) + $newend : int(rand($posbuffer)) + $newend; $newend = $orignewend if $newend <= 0;
		}
		print $out "$chr\t$newstart\t$newend\t$name\t$val\t$strand\t$others\n" if defined($others) and $newstart < $newend;
		print $out "$chr\t$newstart\t$newend\t$name\t$val\t$strand\n" if not defined($others) and $newstart < $newend;
		print $out "$chr\t$newend\t$newstart\t$name\t$val\t$strand\t$others\n" if defined($others) and $newstart >= $newend;
		print $out "$chr\t$newend\t$newstart\t$name\t$val\t$strand\n" if not defined($others) and $newstart >= $newend;
	}
}
close $in;
close $out;
