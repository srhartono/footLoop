#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_m $opt_M $opt_i);
getopts("vM:m:i:");

my ($mem, $max_squeue) = ($opt_m, $opt_M);
$mem = 4000 if not defined $opt_m;
$max_squeue = 50 if not defined $opt_M;
my ($user) = $ENV{HOME} =~ /home\/(\w+)\/$/;
	($user) = $ENV{HOME} =~ /home\/(\w+)/ if not defined $user;
my ($input1) = ($opt_i);
die "\nusage: $YW$0$N $CY-i <input.bam>$N

Options:
-m memory per paralel run [4000] (4000 = 4Gb)
-M maximum squeue to be run [50]

" unless defined $opt_i and -e $opt_i;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");

my $folderFull = getFullpath($folder1);
$folderFull =~ s/\/\.\//\//g;
$folderFull =~ s/\/+/\//g;
print "use=$user, $ENV{HOME}\n";
makedir("$folderFull/sbatch") if not -d "$folderFull/sbatch";
makedir("$folderFull/sbatch/0_CCS2") if not -d "$folderFull/sbatch/0_CCS2";
makedir("$folderFull/sbatch/1_LOG2") if not -d "$folderFull/sbatch/1_LOG2";
my $linecount = 0;
my %data;
my %run;
my @zmw;
my $switch = 0;
my $total_zmw = 0;
my %len;
$len{200000} = 0;
$len{100000} = 0;
$len{50000} = 0;
$len{10000} = 0;
$len{0} = 0;
open (my $in1, "samtools view $input1|") or die "Cannot read from $input1: $!\n";
open (my $out0, ">", "sbatch_record.tsv") or die "Failed to write to sbatch_record.tsv: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 10;
	my ($zmw) = $arr[0] =~ /^.+\/(\d+)\/\d+_\d+$/;
		($zmw) = $arr[0] =~ /^.+\/(\d+)\/\d+_\d+/ if not defined $zmw;
	my $length = length($arr[9]);
	die "can't get zmw from $arr[0]\n" if not defined $zmw;
#	print "$linecount total=$total_zmw ZMW=$zmw\n" if $linecount % 100 == 0;
	if (not defined $data{$zmw}) {
		$total_zmw ++;
		foreach my $currzmw (sort {$a <=> $b} keys %data) {
			my $total_len = $data{$currzmw}; $total_len = 0 if not defined $total_len;
			next if defined $run{$currzmw};
			if ($total_len > 0) {
				$run{$currzmw} = $total_len;
			}
		}
		if ((keys %run) == 100) {
			my $total_lenz = 0;
			my $zmws;
			foreach my $currzmw (sort {$a <=> $b} keys %run) {
				$total_lenz = $run{$currzmw};
				$len{200000} ++ if $total_lenz >= 200000;
				$len{100000} ++ if $total_lenz >= 100000 and $total_lenz <  200000;
				$len{50000} ++ if $total_lenz >= 50000 and $total_lenz <  100000;
				$len{10000} ++ if $total_lenz >= 10000 and $total_lenz <  50000;
				$len{0} ++ if $total_lenz >= 0 and $total_lenz <  10000;
				$zmws .= "$currzmw,";
#				print "$currzmw: $total_lenz\n";
			}
#			foreach my $length (sort {$a <=> $b} keys %len) {
#				print "$length: $len{$length}\n";
#			}
			$zmws =~ s/,$//;
			my $cmd0 = "module load smrtlink/6.0";
			my $cmd  = "time ccs --zmws $zmws --reportFile $folderFull/sbatch/1_LOG2/$fileName1\_$total_zmw\.bam_report.txt --logFile $folderFull/sbatch/1_LOG2/$fileName1\_$total_zmw\.bam_log.txt $folderFull/$fileName1 $folderFull/sbatch/0_CCS2/$fileName1\_$total_zmw\.bam_ccs.fq";
			my $out1;
			open ($out1, ">", "$folderFull/sbatch/sbatch_$total_zmw\.sbatch") or die "Failed to write to $folderFull/sbatch/sbatch_$total_zmw\.sbatch: $!\n";
			print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem $mem -t 99999 --job-name long_sbatch_$total_zmw --output \%j_long_sbatch_$total_zmw.sbout\necho \"total_length=$total_lenz\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
			print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem $mem -t 99999 --job-name long_sbatch_$total_zmw --output \%j_long_sbatch_$total_zmw.sbout\necho \"total_length=$total_lenz\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
			close $out1;
			my ($total_squeue) = `squeue|grep $user|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
			while ($total_squeue > $max_squeue) {
				sleep 60;
				($total_squeue) = `squeue|grep $user|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				print date() . "Sleeping 60s as total squeue is $total_squeue\n";
			}
			my $log = `sbatch $folderFull/sbatch/sbatch_$total_zmw\.sbatch && /bin/rm $folderFull/sbatch/sbatch_$total_zmw\.sbatch`;
			print $out0 ">sbatch_$total_zmw\.sbatch LOG=\n$log\n";
#			print ">sbatch_$total_zmw\.sbatch LOG=\n$log\n";
			$switch = $switch == 0 ? 1 : 0;
			%run = ();
			%data = ();
		}
	}
#	push(@zmw, $zmw) if not grep(/^$zmw$/, @zmw) and not defined $data{$zmw};
	#print "$zmw, total zmw array=" . scalar(@zmw) . "\n" if not defined $data{$zmw};
	$data{$zmw} += $length;
	print "Done $linecount ($total_zmw total zmw, 0=$len{0}, 10k=$len{10000}, 50k=$len{50000}, 100k=$len{100000}, 200k=$len{200000})\n" if $linecount % 100000 == 0;
	$linecount ++;
#	die if $linecount > 50000;
	#next if $total_zmw < 2500;
}
close $in1;


		foreach my $currzmw (sort {$a <=> $b} keys %data) {
			my $total_len = $data{$currzmw};
			next if defined $run{$currzmw};
			if ($total_len > 0) {
				$run{$currzmw} = $total_len;
			}
		}
		if (keys %run != 0) {
			my $total_len = 0;
			my $zmws;
			foreach my $currzmw (sort {$a <=> $b} keys %run) {
				$total_len += $run{$currzmw};
				$zmws .= "$currzmw,";
			}
			$zmws =~ s/,$//;
			my $cmd0 = "module load smrtlink/6.0";
			my $cmd  = "time ccs --zmws $zmws --reportFile $folderFull/sbatch/1_LOG2/$fileName1\_$total_zmw\.bam_report.txt --logFile $folderFull/sbatch/1_LOG2/$fileName1\_$total_zmw\.bam_log.txt $folderFull/$fileName1 $folderFull/sbatch/0_CCS2/$fileName1\_$total_zmw\.bam_ccs.fq";
			my $out1;
			open ($out1, ">", "$folderFull/sbatch/sbatch_$total_zmw\.sbatch") or die "Failed to write to $folderFull/sbatch/sbatch_$total_zmw\.sbatch: $!\n";
			print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem $mem -t 99999 --job-name long_sbatch_$total_zmw --output \%j_long_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
			print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem $mem -t 99999 --job-name long_sbatch_$total_zmw --output \%j_long_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
			close $out1;
			my ($total_squeue) = `squeue|grep $user|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
			while ($total_squeue > $max_squeue) {
				sleep 60;
				($total_squeue) = `squeue|grep $user|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				print date() . "Sleeping 60s as total squeue is $total_squeue\n";
			}
			my $log = `sbatch $folderFull/sbatch/sbatch_$total_zmw\.sbatch && /bin/rm $folderFull/sbatch/sbatch_$total_zmw\.sbatch`;
			print $out0 ">sbatch_$total_zmw\.sbatch LOG=\n$log\n";
			$switch = $switch == 0 ? 1 : 0;
			%run = ();
			%data = ();
		}
