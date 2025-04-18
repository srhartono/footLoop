#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use colorz;
use vars qw($opt_n $opt_0);
getopts("n:0");


my ($footPeakFolder) = $opt_n;
($footPeakFolder) = @ARGV if not defined $opt_n;
die "\nUsage: $YW$0$N <dir> or -n <footPeakfolder>\n\n" if not defined $footPeakFolder;
die "\nUsage: $YW$0$N -n <footPeakfolder>\nfootPeakFolder doesn't exist! ($LCY$footPeakFolder$N)\n\n" if defined $footPeakFolder and not -e $footPeakFolder;
my @files;
my $goodFolder = "goodFolderUNK";
if (@files == 0) {
	@files = <$footPeakFolder/.CALL/*.[NP][OE][PA]K.out>;
	$goodFolder = "$footPeakFolder/.CALL/";
	print "\nfiles = $LCY<$footPeakFolder/.CALL/*.[NP][OE][PA]K.out>$N\n\n" if @files > 0;
}
if (@files == 0) {
	@files = <$footPeakFolder/*.[NP][OE][PA]K.out>;
	print "\nfiles = $LCY<$footPeakFolder/*.[NP][OE][PA]K.out>$N\n\n" if @files > 0;
	$goodFolder = "$footPeakFolder/";
}
my $totalfile = @files;
if (@files == 0) {
	print "There is ${LCY}0$N *.out files in $footPeakFolder/ or $footPeakFolder/.CALL/\n\n";
	exit 0;
}
#@files = <./*.[NP][OE][PA]K.out> if @files == 0;
print "Processing $LCY$goodFolder$N\n";
my @cmd;
my $cmdcount = 0;
my $todofile = 0;
for (my $i = 0; $i < @files; $i++) {
	my $iprint = $i+1;
	my $file = $files[$i];
	my $wclfile = $file;
	$wclfile =~ s/.out$/.wcl/;
	print "$YW$iprint/$totalfile$N $cmdcount $LGN$wclfile$N\n" if $i == 0 or ($i+1) % 1000 == 0 or $i == $totalfile - 1;
	next if -e $wclfile and -s $wclfile > 0;
	$todofile ++;
	print "$YW$iprint/$totalfile$N $cmdcount Todo #$LGN$todofile$N $LCY$wclfile$N\n" if $todofile % 1000 == 0;
	my $cmd = "wc -l $file > $wclfile";
	if (@cmd < 250) {
		push(@cmd, $cmd);
	}
	else {
		my $outwclsbatch = "wclthese.$cmdcount.sbatch";
		$cmdcount ++;
		print "$outwclsbatch\n";
		if (not defined $opt_0) {
		open (my $out1, ">", "wclthese.$cmdcount.sbatch") or die;
		print $out1 "#!/bin/bash -l
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p high
#SBATCH --mem 2000
#SBATCH --job-name \"wclthese $cmdcount\"
#SBATCH --output \".%j_wclthese.$cmdcount\"
#SBATCH -t 999:99:99\n
";
		print $out1 join("\n", @cmd);
		close $out1;
		system("sbatch $outwclsbatch") == 0 or die "\n\nFailed to sbatch $LCY$outwclsbatch$N: $!\n\n" if not defined $opt_0;
		}
		@cmd = ();
	}
}
if (@cmd > 0 and not defined $opt_0) {
	my $outwclsbatch = "wclthese.$cmdcount.sbatch";
	$cmdcount ++;
	if (not defined $opt_0) {
	open (my $out1, ">", "wclthese.$cmdcount.sbatch") or die;
		print $out1 "#!/bin/bash -l
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p high
#SBATCH --mem 2000
#SBATCH --job-name \"wclthese $cmdcount\"
#SBATCH --output \".%j_wclthese.$cmdcount\"
#SBATCH -t 999:99:99\n
";
	print $out1 join("\n", @cmd);
	close $out1;
	system("sbatch $outwclsbatch") == 0 or die "\n\nFailed to sbatch $LCY$outwclsbatch$N: $!\n\n" if not defined $opt_0;
	}
	@cmd = ();
}
print "\n\ntotal file = $totalfile, todofile=$todofile\n\n";
print "There is $LGN$totalfile$N *.out files in $LGN$goodFolder$N ($LGN$cmdcount$N sbatch files)\n\n";
