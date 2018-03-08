#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use Thread; use Thread::Queue;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_n);
getopts("vn:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';

my $date = getDate();
my $uuid = getuuid();

die "\nusage: $YW$0$N -n $CY<FootPeakFolder>$N\n\n" unless defined $opt_n and -d $opt_n;
my ($footPeakFolder) = (getFullpath($opt_n));
my @files = (<$footPeakFolder/.CALL/*.PEAK>, <$footPeakFolder/.CALL/*.NOPK>);

print "footPeakFolder = $footPeakFolder\n";
foreach my $file (sort @files) {
	my ($folder1, $fileName1) = mitochy::getFilename($file, "folderfull");
	my ($gene, $strand, $window, $thres, $type) = $fileName1 =~ /^(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC)/;
	my $peakFile = $folder . "/$gene\_$strand\_$window\_$thres\_$type.PEAK";
	
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

