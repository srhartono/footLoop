#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v);
getopts("v");

my ($input1) = @ARGV;
die "\nusage: $YW$0$N $CY<input1>$N\n\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my $user = 0; #0 mitochi, 1 shanaya
if ($ENV{HOME} =~ /mitochi/) {
	$user = 1;
}
print "use=$user, $ENV{HOME}\n";
if ($user eq 0) {
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/";
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/0_CCS2") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/0_CCS2";
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/1_LOG2") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/1_LOG2";
}
else {
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/";
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/0_CCS2") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/0_CCS2";
	makedir("/home/mitochi/shared/shah/seq1811/sbatch/1_LOG2") if not -d "/home/mitochi/shared/shah/seq1811/sbatch/1_LOG2";
}
my $linecount = 0;
my %data;
my %run;
my @zmw;
my $switch = 0;
my $total_zmw;
my %len;
$len{200000} = 0;
$len{100000} = 0;
$len{50000} = 0;
$len{10000} = 0;
$len{0} = 0;
open (my $in1, "samtools view $input1|") or die "Cannot read from $input1: $!\n";
open (my $out0, ">", "sbatch_record_longboi_evennumber.tsv") or die "Failed to write to sbatch_record_longboi_evennumber.tsv: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 10;
	my ($zmw) = $arr[0] =~ /^.+\/(\d+)\/\d+_\d+$/;
	my $length = length($arr[9]);
	die "can't get zmw from $arr[0]\n" if not defined $zmw;
	if (not defined $data{$zmw}) {
		$total_zmw ++;
		foreach my $currzmw (sort {$a <=> $b} keys %data) {
			my $total_len = $data{$currzmw}; $total_len = 0 if not defined $total_len;
			next if defined $run{$currzmw};
			if ($total_len >= 50000) {
				$run{$currzmw} = $total_len;
			}
		}
		if (keys %run == 1000) {
			my $total_len = 0;
			my $zmws;
			foreach my $currzmw (sort {$a <=> $b} keys %run) {
					$total_len += $run{$currzmw};
					$len{200000} ++ if $total_len >= 200000;
					$len{100000} ++ if $total_len >= 100000 and $total_len <  200000;
					$len{50000} ++ if $total_len >= 50000 and $total_len <  100000;
					$len{10000} ++ if $total_len >= 10000 and $total_len <  50000;
					$len{0} ++ if $total_len >= 0 and $total_len <  10000;
					$zmws .= "$currzmw,";
			}
			$zmws =~ s/,$//;
			my $cmd0 = "module load smrtlink/6.0";
			my $cmd  = "time ccs --zmws $zmws --reportFile /home/mitochi/shared/shah/seq1811/sbatch/1_LOG2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_report.txt --logFile /home/mitochi/shared/shah/seq1811/sbatch/1_LOG2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_log.txt /home/mitochi/shared/shah/seq1811/m54290_181111_201646.subreads.bam /home/mitochi/shared/shah/seq1811/sbatch/0_CCS2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_ccs.fq";
			my $out1;
			if ($user == 0 and $switch == 0) {
				die;
				open ($out1, ">", "/home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch") or die "Failed to write to /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch: $!\n";
				print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_oddnumber_sbatch_$total_zmw --output \%j_long_oddnumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_oddnumber_sbatch_$total_zmw --output \%j_long_oddnumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				close $out1;
				my ($total_squeue) = `squeue|grep mitochi|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				while ($total_squeue > 50) {
					sleep 60;
					($total_squeue) = `squeue|grep mitochi|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				}
				my $log = `sbatch /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch && /bin/rm /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch`;
				print $out0 ">sbatch_$total_zmw\_longboi_oddnumber.sbatch LOG=\n$log\n";
			}
			elsif ($user == 1 and $switch == 1) {
				$cmd =~ s/mitochi/mitochi/g;
				open ($out1, ">", "/home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch") or die "Failed to write to /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch: $!\n";
				print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_evennumber_sbatch_$total_zmw --output \%j_long_evennumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_evennumber_sbatch_$total_zmw --output \%j_long_evennumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				close $out1;
				my ($total_squeue) = `squeue|grep mitochi|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				while ($total_squeue > 50) {
					sleep 60;
						($total_squeue) = `squeue|grep mitochi|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				}
				my $log = `sbatch /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch && /bin/rm /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch`;
				print $out0 ">sbatch_$total_zmw\_longboi_evennumber.sbatch LOG=\n$log\n";
			}
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
	#next if $total_zmw < 2500;
}
close $in1;


		foreach my $currzmw (sort {$a <=> $b} keys %data) {
			my $total_len = $data{$currzmw};
			next if defined $run{$currzmw};
			if ($total_len >= 50000) {
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
			my $cmd  = "time ccs --zmws $zmws --reportFile /home/mitochi/shared/shah/seq1811/sbatch/1_LOG2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_report.txt --logFile /home/mitochi/shared/shah/seq1811/sbatch/1_LOG2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_log.txt /home/mitochi/shared/shah/seq1811/m54290_181111_201646.subreads.bam /home/mitochi/shared/shah/seq1811/sbatch/0_CCS2/m54290_181111_201646.subreads_$total_zmw\_longboi_evennumber.bam_ccs.fq";
			my $out1;
			if ($user == 0 and $switch == 0) {
				open ($out1, ">", "/home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch") or die "Failed to write to /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch: $!\n";
				print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_oddnumber_sbatch_$total_zmw --output \%j_long_oddnumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_oddnumber_sbatch_$total_zmw --output \%j_long_oddnumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				close $out1;
				my ($total_squeue) = `squeue|grep mitochi|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				while ($total_squeue > 1000) {
					sleep 60;
					($total_squeue) = `squeue|grep mitochi|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				}
				my $log = `sbatch /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch && /bin/rm /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_oddnumber.sbatch`;
				print $out0 ">sbatch_$total_zmw\_longboi_oddnumber.sbatch LOG=\n$log\n";
			}
			elsif ($user == 1 and $switch == 1) {
				$cmd =~ s/mitochi/mitochi/g;
				open ($out1, ">", "/home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch") or die "Failed to write to /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch: $!\n";
				print $out0 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_evennumber_sbatch_$total_zmw --output \%j_long_evennumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				print $out1 "#!/bin/bash -l\n#SBATCH -p high --mem 10000 -t 99999 --job-name long_evennumber_sbatch_$total_zmw --output \%j_long_evennumber_sbatch_$total_zmw.sbout\necho \"total_length=$total_len\"\necho \"$cmd0\"\necho \"$cmd\"\n$cmd0\n$cmd\n\n";
				close $out1;
				my ($total_squeue) = `squeue|grep mitochi|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				while ($total_squeue > 1000) {
					sleep 60;
						($total_squeue) = `squeue|grep mitochi|grep -v bigmem|wc -l` =~ /^(\d+)/; $total_squeue = 0 if not defined $total_squeue; chomp($total_squeue);
				}
				my $log = `sbatch /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch && /bin/rm /home/mitochi/shared/shah/seq1811/sbatch/sbatch_$total_zmw\_longboi_evennumber.sbatch`;
				print $out0 ">sbatch_$total_zmw\_longboi_evennumber.sbatch LOG=\n$log\n";
			}
			$switch = $switch == 0 ? 1 : 0;
			%run = ();
			%data = ();
		}
