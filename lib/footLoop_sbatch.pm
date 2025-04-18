package footLoop_sbatch;

# footLoop pipeline Version 3
# Copyright (C) 2019-2029 Stella Regina Hartono
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The license can be found at https://www.gnu.org/licenses/gpl-3.0.en.html.
# By downloading or using this software, you agree to the terms and conditions of the license.

use strict; use warnings;
use vars qw(@EXPORT);
use parent "Exporter";

our @EXPORT = qw(footLoop_sbatch_main squeue_check print_cmd);

use myFootLib;
use FAlite;

#my $homedir = $ENV{"HOME"};
#my ($footLoop_script_folder, $version, $md5script) = check_software();
#my $footLoopScriptsFolder = dirname(dirname abs_path $0);
#my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
#$version_small = $version if not defined $version_small;

sub footLoop_sbatch_main {
	my ($cmdtemplate, $suffix, $filesARRAY, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir, $mem, $debug) = @_;
	my $force_sbatch_print = defined $force_sbatch ? "force_sbatch is on" : "force_sbatch_is_off";
	LOG($outLog, "${LPR}footLoop_sbatch::footLoop_sbatch_main$N $force_sbatch_print\n");
	if (defined $debug) {
		LOG($outLog, "${LPR}footLoop_sbatch::footLoop_sbatch_main$N debug is on!\n");
	}
	else {
		LOG($outLog, "${LPR}footLoop_sbatch::footLoop_sbatch_main$N debug is off!\n");
	}
	
	my %force;
	
	# - suffix : $file\_$suffix.sbatch/sbout/folder
	# - ext	 : determine donefile (<$file\_$suffix/*.$ext>)

	my $jobidhash;
	my $totalfiles = scalar(@{$filesARRAY});
	for (my $i = 0; $i < @{$filesARRAY}; $i++) {
		my $iprint = $i + 1;
		my $file = $filesARRAY->[$i];
		my ($folder, $filename) = getFilename($file, "folderfull");
		$outsbatchDir = $folder if not defined $outsbatchDir;
		my $sbatchfile = "$outsbatchDir/$filename\_$suffix.sbatch";
		my $sboutfile  = "$outsbatchDir/$filename\_$suffix.sbout";
		my $donefile	= "$outsbatchDir/$filename\_$suffix.done";
		
		my $cmd = $cmdtemplate;
			$cmd =~ s/FILENAME/$file/g;	

		if ($cmd !~ /FOLDER/) {
			system("mkdir -p $outsbatchDir/$filename\_$suffix") if not -e "$outsbatchDir/$filename\_$suffix";
		}
		else {
			$cmd =~ s/FOLDER/$outsbatchDir/g;	
		}
		if ($i == 0) {
			print_cmd($cmd, $outLog);
		}

		if ($cmdtemplate =~ /FNINDICE/) {
			$cmd =~ s/FNINDICE/$i/g;
		}

		$sboutfile =~ s/\/+/\//g;
		$mem = 16000 if not defined $mem;
		my $sbatchprint = "";
			$sbatchprint .= "#!/bin/bash -l\n";
			$sbatchprint .= "#SBATCH -n 2 -N 1 -p high --mem $mem -t 999:99:99\n";
			$sbatchprint .= "#SBATCH --job-name \"$filename\_$suffix\"\n";
			$sbatchprint .= "#SBATCH --output \"$sboutfile\"\n\n";
			$sbatchprint .= "conda activate footLoop2\n";
			$sbatchprint .= "$cmd && echo \"Done!\" > $donefile\n\n";

		#"sbatchFile=\n$LCY$sbatchfile$N\n\n" if defined $debug;
		open (my $out, ">", $sbatchfile) or die "Can't write to $LCY$sbatchfile$N: $!\n";
		print $out $sbatchprint;
		close $out;
		my $is_run =  (not defined $debug) ? "$LGN(PRINTED AND RAN)$N" : "$LPR(PRINTED BUT NOT RUN/DEBUG)$N";
		if (not defined $force{0} and not defined $force_sbatch and -e $donefile) {
			LOG($outLog, "\n" . date() . "${LPR}$iprint/$totalfiles sbatch_these $suffix$N: sbatch $LCY$sbatchfile$N # ${LGN}DONE$N\n");
			LOG($outLog, date() . "${LPR}$iprint/$totalfiles donefile: $LCY$donefile$N\n");
			next;
		}
		else {
			LOG($outLog, "\n" . date() . "${LPR}$iprint/$totalfiles sbatch_these $suffix$N: sbatch: $LCY$sbatchfile$N $is_run\n");
			LOG($outLog, date() . "${LPR}$iprint/$totalfiles donefile: $LCY$donefile$N\n");
			#LOG($outLog, date() . "${LPR}$iprint/$totalfiles donefile: $LCY$donefile$N\n") if defined $debug;
		}
		
		if (defined $debug) { # Debug
			next;
		}

		if ($i != 0) {
			my $sleep = 0;
			while (1) {
				last if $i < $max_parallel_run;
				my ($job_left, $test_sboutfile) = squeue_check($jobidhash);
				LOG($outLog, "\n" . date() . "[A] $job_left jobs left! ($LCY$test_sboutfile$N)\n") if $sleep % 12 == 0;
				last if ($job_left < $max_parallel_run);
				$sleep ++;
				sleep 5;
			}
		}
		my ($jobid) = `sbatch $sbatchfile`;
		chomp($jobid);
		($jobid) = $jobid =~ /^Submi.+job (\d+)$/;
		next if not defined $jobid;
		$jobidhash->{$jobid} = $sboutfile;

		LOG($outLog, "$YW$i$N $LCY$filename$N $LGN$jobid$N\n");
	}
	my $sleep = 0;
	while (1) {
		my ($job_left, $test_sboutfile) = squeue_check($jobidhash);
		my $sleeptime = $job_left < $max_parallel_run ? 5 : $job_left < 100 ? 20 : 60;
		LOG($outLog, "\n" . date() . "[B] $job_left jobs left! ($LCY$test_sboutfile$N)\n") if $sleep % 12 == 0;
		#LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
		last if $job_left == 0;
		$sleep ++;
		sleep $sleeptime;
	}
	LOG($outLog, "\n" . date() . "All have been run!\n\n");
	return(0);
}

sub print_cmd {
	my ($cmd, $outLog) = @_;
	if ($cmd !~ /^#/) {
		$cmd =~ s/^/    /;
		$cmd =~ s/ \-/ \\\n      \-/g;
	}
	else {
		$cmd =~ s/^#/   # /;
		$cmd =~ s/ \-/ \\\n     # \-/g;
	}
	LOG($outLog, "$LGN\n$cmd\n$N\n");
}

sub squeue_check {
	my ($jobidhash, $outLog) = @_;
	my $test_sboutfile;
	my @squeue = `squeue`;
	my $squeuehash;
	foreach my $line (@squeue) {
		next if $line =~ /JOBID\s+PARTITION.+/;
		my ($jobid) = $line =~ /^\s*(\d+)\s+/;
		if (not defined $jobid) {
			LOG($outLog, "Can't parse jobid from line=$LCY$line$N\n");
			next; # just next so we don't kill the script...
		}
		next if not defined $jobidhash->{$jobid};
		$squeuehash->{$jobid} = 1;
	}
	foreach my $jobid (keys %{$jobidhash}) {
		if (defined $squeuehash->{$jobid}) {
			$test_sboutfile = $jobidhash->{$jobid} if not defined $test_sboutfile;
			next;
		}
		undef $jobidhash->{$jobid};
		delete $jobidhash->{$jobid};
	}
	my ($total) = scalar(keys %{$jobidhash});
	$test_sboutfile = "[none]" if not defined $test_sboutfile;
	return ($total, $test_sboutfile);
}
