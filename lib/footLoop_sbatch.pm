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

our @EXPORT = qw(footLoop_sbatch_main);

use myFootLib;
use FAlite;

#my $homedir = $ENV{"HOME"};
#my ($footLoop_script_folder, $version, $md5script) = check_software();
#my $footLoopScriptsFolder = dirname(dirname abs_path $0);
#my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
#$version_small = $version if not defined $version_small;

sub footLoop_sbatch_main {
	my ($cmdtemplate, $suffix, $filesARRAY, $max_parallel_run, $outLog, $force_sbatch, $folderwant, $debug) = @_;

	my %force;
	
	# - suffix : $file\_$suffix.sbatch/sbout/folder
	# - ext	 : determine donefile (<$file\_$suffix/*.$ext>)

	my $jobidhash;
	my $totalfiles = scalar(@{$filesARRAY});
	for (my $i = 0; $i < @{$filesARRAY}; $i++) {

		my $file = $filesARRAY->[$i];
		my ($folder, $filename) = getFilename($file, "folderfull");
		$folderwant = $folder if not defined $folderwant;
		my $sbatchfile = "$folderwant/$filename\_$suffix.sbatch";
		my $sboutfile  = "$folderwant/$filename\_$suffix.sbout";
		my $donefile	= "$folderwant/$filename\_$suffix.done";
		
		my $cmd = $cmdtemplate;
			$cmd =~ s/FILENAME/$file/g;	

		if ($cmd !~ /FOLDER/) {
			system("mkdir -p $folderwant/$filename\_$suffix") if not -e "$folderwant/$filename\_$suffix";
		}
		else {
			$cmd =~ s/FOLDER/$folderwant/g;	
		}
		if ($i == 0) {
			print_cmd($cmd, $outLog);
		}

		if ($cmdtemplate =~ /FNINDICE/) {
			$cmd =~ s/FNINDICE/$i/g;
		}

		$sboutfile =~ s/\/+/\//g;
		my $sbatchprint = "";
			$sbatchprint .= "#!/bin/bash -l\n";
			$sbatchprint .= "#SBATCH -n 2 -N 1 -p high --mem 16000 -t 999:99:99\n";
			$sbatchprint .= "#SBATCH --job-name \"$filename\_$suffix\"\n";
			$sbatchprint .= "#SBATCH --output \"$sboutfile\"\n\n";
			$sbatchprint .= "conda activate footLoop2\n";
			$sbatchprint .= "$cmd && echo \"Done!\" > $donefile\n\n";

		#"sbatchFile=\n$LCY$sbatchfile$N\n\n" if defined $debug;
		open (my $out, ">", $sbatchfile) or die "Can't write to $LCY$sbatchfile$N: $!\n";
		print $out $sbatchprint;
		close $out;

		if (not defined $force{0} and not defined $force_sbatch and -e $donefile) {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles sbatch_these $suffix$N: sbatch $LCY$sbatchfile$N # ${LGN}DONE$N\n");
			next;
		}
		else {
			LOG($outLog, "\n" . date() . "${LPR}$i/$totalfiles sbatch_these $suffix$N: sbatch: $LCY$sbatchfile$N\n");
		}
		
		if (defined $debug) { # Debug
			next;
		}

		if ($i != 0) {
			my $sleep = 0;
			while (1) {
				last if $i < $max_parallel_run;
				my ($job_left) = squeue_check($jobidhash);
				LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
				last if ($job_left < $max_parallel_run);
				$sleep ++;
				sleep 5;
			}
		}
		my ($jobid) = `sbatch $sbatchfile`;
		chomp($jobid);
		($jobid) = $jobid =~ /^Submi.+job (\d+)$/;
		next if not defined $jobid;
		$jobidhash->{$jobid} = 1;
		LOG($outLog, "$YW$i$N $LCY$filename$N $LGN$jobid$N\n");
	}
	my $sleep = 0;
	while (1) {
		my ($job_left) = squeue_check($jobidhash);
		my $sleeptime = $job_left < $max_parallel_run ? 5 : $job_left < 100 ? 20 : 60;
		LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
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
		next if defined $squeuehash->{$jobid};
		undef $jobidhash->{$jobid};
		delete $jobidhash->{$jobid};
	}
	my ($total) = scalar(keys %{$jobidhash});
	return ($total);
}
