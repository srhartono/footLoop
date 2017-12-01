#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_x $opt_R);
getopts("vxR");

my ($folders, $lotsOfC) = @ARGV;
die "\nusage: $YW$0$N $CY<input1_Pos50.txt>$N\n\n" unless @ARGV;
my %bad;
if (defined $lotsOfC) {
	my @lotsOfC = split(";", $lotsOfC);
	foreach my $coor (@lotsOfC) {
		my ($beg, $end) = split(",", $coor);
		$bad{$beg} = $end;
	}
}
my %pk;
my (@inputs, $input1);
if ($folders =~ /.txt$/) {
	$input1 = $folders;
}
else {
	(@inputs) = <$folders/*Pos50.txt>;
	die "Must have 1 input only! (" . scalar(@inputs) . "):\n" . join("\n-", @inputs) . "\n" if @inputs != 1 and not defined $opt_x;
	$input1 = defined $opt_x ? $folders : $inputs[0];
}
my ($folder, $fileName) = getFilename($input1, "folder");
my %total; 
$total{Pos}{peak} = 0; $total{Pos}{nopeak} = 0; $total{Pos}{total} = 0;
$total{Neg}{peak} = 0; $total{Neg}{nopeak} = 0; $total{Neg}{total} = 0;
my $type = "Pos";
my $log2 = "";
print STDERR "Folder $YW$folder$N: Processing files related to $LCY$input1$N\n";
open (my $outLog, ">", "$folder/0_RESULTS.TXT") if not defined $opt_x;
open (my $outLog2, ">", "$folder/1_RESULTS_EXTRA.TXT") if not defined $opt_x;
for (my $h = 0; $h < 2; $h++) {
	my $peakFile   = $input1; $peakFile   =~ s/_NOPEAK.txt/.txt/ if $peakFile =~ /_NOPEAK.txt/;
	my $nopeakFile = $input1; $nopeakFile =~ s/.txt$/_NOPEAK.txt/ if $nopeakFile !~ /_NOPEAK.txt/;
	if (not defined $opt_x) {
		$peakFile =~ s/_Neg/_Pos/ if $h == 0;
		$peakFile =~ s/_Pos/_Neg/ if $h == 1;
		$nopeakFile =~ s/_Neg/_Pos/ if $h == 0;
		$nopeakFile =~ s/_Pos/_Neg/ if $h == 1;
	}
	last if $h == 1 and defined $opt_x;
	$type = $h == 0 ? "Pos" : "Neg";
	print STDERR "\t$h\t$peakFile\t$nopeakFile\n";
#	die if $h == 1;
#	next;
	
	my ($folder1, $peakfileName) = mitochy::getFilename($peakFile, "folder");
	my ($folder2, $nopeakfileName) = mitochy::getFilename($nopeakFile, "folder");
	
	my %data; my $totalnopeak = 0;
	open (my $in2, "<", $nopeakFile) or die "Cannot read from $nopeakFile: $!\n";
	my ($totalline) = `wc -l $nopeakFile` =~ /^(\d+)/;
	print STDERR "\tProcessing $LPR$nopeakFile$N\n";
	while (my $line = <$in2>) {
		chomp($line);
		print STDERR "\t$CY$nopeakFile$N: Done $totalnopeak / $totalline\n" if $totalnopeak % 500 == 0;
		next if $line =~ /^#/;
		my ($name, $val, $totalPeak, $peaks) = parse_peak($line);
		$val = "$name\t" . join("\t", @{$val});
		push(@{$data{peak}}, $val) if $totalPeak > 0;
		push(@{$data{nopeak}}, $val) if $totalPeak == 0;
		$totalnopeak ++;
		$pk{$peakFile}{$name} = $peaks if defined $peaks;
	}
	close $in2;
	
	my $peakCount = defined $data{peak} ? @{$data{peak}} : 0;
	my $nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} : 0;
	my $nopeakPrint ="$folder2\t$nopeakfileName\t$peakCount\t$nopeakCount\t$totalnopeak\t$totalline";
	$total{$type}{peak} += $peakCount;
	$total{$type}{nopeak} += $nopeakCount;
	$total{$type}{total} += $totalnopeak;
	
	my $totalpeak = 0;
	open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";
	($totalline) = `wc -l $peakFile` =~ /^(\d+)/;
	print STDERR "\tProcessing $LCY$peakFile$N\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		print STDERR "\t$CY$peakFile$N: Done $totalpeak / $totalline\n" if $totalpeak % 100 == 0;
		my ($name, $val, $totalPeak, $peaks) = parse_peak($line);
		$val = "$name\t" . join("\t", @{$val});
		push(@{$data{peak}}, $val) if $totalPeak > 0;
		push(@{$data{nopeak}}, $val) if $totalPeak == 0;
		$totalpeak ++;
		$pk{$peakFile}{$name} = $peaks if defined $peaks;
	}
	close $in1;
	
	$peakCount = defined $data{peak} ? @{$data{peak}} - $peakCount : 0;
	$nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} - $nopeakCount : 0;
	my $peakPrint ="$folder1\t$peakfileName\t$peakCount\t$nopeakCount\t$totalpeak\t$totalline";
	$total{$type}{peak} += $peakCount;
	$total{$type}{nopeak} += $nopeakCount;
	$total{$type}{total} += $totalpeak;



	open (my $out1, ">", "$folder1/$peakfileName.out") or die "Cannot write to $peakfileName.out: $!\n";
	open (my $out2, ">", "$folder1/$nopeakfileName.out") or die "Cannot write to $nopeakfileName.out: $!\n";
	foreach my $val (sort @{$data{peak}}) {
		print $out1 "$val\n";
	}
	foreach my $val (sort @{$data{nopeak}}) {
		print $out2 "$val\n";
	}
	
	close $out1;
	close $out2;
#### HERE ##
#	next;
	
	$log2 .= "\#Folder\tFile\tPeak\tNoPeak\tTotalRead\tTotalLineInFile\n" if $h == 0 and not defined $opt_x;
	print STDERR "$peakPrint\n$nopeakPrint\n" if defined $opt_x;
	$log2 .= "$peakPrint\n$nopeakPrint\n" if not defined $opt_x;

	mkdir "$folder1/remove" if not -d "$folder1/remove/";
	
	my $peakFileBackup = "$folder1/remove/$peakfileName.txt";
	my $peakFileTemp   = $peakFileBackup;
	my $count = 0;
	while (-e $peakFileTemp) {
		$peakFileTemp = $peakFileBackup . $count;
		$count ++;
	}
	print STDERR "\tmv $peakFile $peakFileTemp\n";
	system("/bin/mv $peakFile $peakFileTemp") if not defined $opt_x;
	print STDERR "\tmv $folder1/$peakfileName.out $folder1/$peakfileName.txt\n";
	system("mv $folder1/$peakfileName.out $folder1/$peakfileName.txt") if not defined $opt_x;
	
	my $nopeakFileBackup = "$folder1/remove/$nopeakfileName.txt";
	my $nopeakFileTemp   = $nopeakFileBackup;
	$count = 0;
	while (-e $nopeakFileTemp) {
		$nopeakFileTemp = $nopeakFileBackup . $count;
		$count ++;
	}
	print STDERR "\t/bin/mv $nopeakFile $nopeakFileTemp\n";
	system("/bin/mv $nopeakFile $nopeakFileTemp") if not defined $opt_x;
	print STDERR "\tmv $folder1/$nopeakfileName.out $folder1/$nopeakfileName.txt\n";
	system("mv $folder1/$nopeakfileName.out $folder1/$nopeakfileName.txt") if not defined $opt_x;
}
#		foreach my $peak (sort @{$peak{peak}}) {
#			my ($beg, $end) = split("-", $peak);

foreach my $file (sort keys %pk) {
	next if not defined $pk{$file};
	next if defined $pk{$file} and keys %{$pk{$file}} == 0;
	open (my $outR, ">", "$file.PEAKS") or die;
	foreach my $name (sort keys %{$pk{$file}}) {
		next if not defined $pk{$file}{$name};
		my @arr = @{$pk{$file}{$name}};
		foreach my $peakz (sort @{$pk{$file}{$name}}) {
			my ($beg, $end) = split("-", $peakz);
			my $name2 = $name; $name2 =~ s/^SEQ_//;
			print $outR "$name2\t$beg\t$end\n";
		}
	}
	close $outR;
}

if (not defined $opt_x) {
	$total{Pos}{total} = 1 if $total{Pos}{total} == 0;
	$total{Pos}{peak} = int(1000 * $total{Pos}{peak} / $total{Pos}{total}+0.5)/10;
	$total{Pos}{nopeak} = int(1000 * $total{Pos}{nopeak} / $total{Pos}{total}+0.5)/10;
	$total{Neg}{total} = 1 if $total{Neg}{total} == 0;
	$total{Neg}{peak} = int(1000 * $total{Neg}{peak} / $total{Neg}{total}+0.5)/10;
	$total{Neg}{nopeak} = int(1000 * $total{Neg}{nopeak} / $total{Neg}{total}+0.5)/10;
	my @folder = split("/", $folder);
	my $foldershort = $folder[@folder-1];
	   $foldershort = $folder[@folder-2] if not defined ($foldershort) or (defined $foldershort and $foldershort =~ /^[\s]*$/);
	my $fileNamePos = $fileName; $fileName =~ s/Neg/Pos/g if $fileName =~ /Neg/;
	my $fileNameNeg = $fileName; $fileName =~ s/Pos/Neg/g if $fileName =~ /Pos/;
	print $outLog "#folder\tfilename\tStrand\ttotal\tpeak.perc\n";
	print $outLog "$foldershort\t$fileNamePos\t$total{Pos}{total}\t$total{Pos}{peak}\n";
	print $outLog "$foldershort\t$fileNameNeg\t$total{Neg}{total}\t$total{Neg}{peak}\n";
	print $outLog2 "$log2";
	close $outLog;
	close $outLog2;
}
system("cat $folder/0_RESULTS.TXT") if not defined $opt_x;
print STDERR "\tcd $folder && run_Rscript.pl *MakeHeatmap.R\n";
system("cd $folder && run_Rscript.pl *MakeHeatmap.R") if not defined $opt_x and defined $opt_R;

###############
# Subroutines #
###############

sub parse_peak {
	my ($name, @val) = split("\t", $_[0]);
	my $peaks;
	my %peak; $peak{curr} = 0; my $edge = 0; my $edge2 = 0; my $six = 0; my $edge1 = 0;
	my $Length = @val; 
	my $print = "Total length = $Length\n";
	for (my $i = 0; $i < @val; $i++) {
		if ($i % 100 == 0 and $edge2 <= 1) {$print .= "\n$YW" . $i . "$N:\t";}
		my $val = $val[$i];
		$six ++ if $val == 6;
		$edge = 1 if $six > 20;
		$six = 0 if $val != 6;
		if ($edge == 1 and $six == 0) {
			$edge = $i;
			$edge1 = $i;
			$print .= "EDGE1 = $edge\n";
			$edge2  = 1;
		}
		if ($edge2 == 1 and $six >= 20) {
			$edge2 = $i - $six;
			$print .= "EDGE2 = $edge2\n";
		}
		if ($val =~ /^(5|9)$/) { # Peak Converted CpG or CH
			$peak{curr} = 1;
			$print .= "${LRD}$val$N";
			if (not defined $peak{beg}) {
				$peak{beg} = $i;
				$peak{end} = $i;
			}
			elsif (defined $peak{beg} and $i - $peak{end} >= 250) {
				push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
				$peak{beg} = $i;
				$peak{end} = $i;
			}
			elsif (defined $peak{beg} and $i - $peak{end} < 250) {
				$peak{end} = $i;
			}
		}
		else {
			$print .= "${LGN}$val$N" if $val =~ /^(1|4)$/;
			$print .= "${LGN}$val$N" if $val =~ /^(0|3)$/;
			$print .= "." if $val =~ /^2$/;
			$print .= "x" if $val =~ /^6$/ and ($edge2 <= 1 or ($edge2 > 1 and $i < $edge2));
			$print .= "EDGE2\n" if $val =~ /^6$/ and ($edge2 > 1 and $i == $edge2);
			if ($peak{curr} == 1 and $i - $peak{end} >= 250) {
				push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
				undef $peak{beg}; undef $peak{end};
				$peak{curr} = 0;
			}
		}
	}
	if ($peak{curr} == 1 and defined $peak{beg}) {
		push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
	}
	my (%nopeak, @peak);
	my %peak2;
	$print .= "\n";
#	print "\nDoing $YW$name$N\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
	if (defined $peak{peak}) {
		foreach my $peak (sort @{$peak{peak}}) {
			my ($beg, $end) = split("-", $peak);
	###		print "$name: $beg to $end\n" if $name eq 77011 or $name eq "SEQ_77011";
			my $checkBad = 0;
			foreach my $begBad (sort keys %bad) {
				my $endBad = $bad{$begBad};
###				print "\t$beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "SEQ_77011";
				if ($beg >= $begBad and $beg <= $endBad and $end >= $begBad and $end <= $endBad) {
#					for (my $m = $beg; $m <= $end; $m++) {
#						$nopeak{$m} = 1;
#					}
###					print "\t\t$LGN YES$N\n" if $name eq "SEQ_77011";
					$checkBad = 1; last;
				}
			}
			if ($checkBad != 1) {
				foreach my $begBad (sort keys %bad) {
					my $endBad = $bad{$begBad};
#					print "$name: $beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
					#if (($beg >= $begBad and $beg <= $endBad) or ($end >= $begBad and $end <= $endBad)) {
						my @valz = @val;
						my ($goodC, $badC) = (0,0);
						for (my $m = $beg; $m < $end; $m++) {
							if ($m >= $begBad and $m <= $endBad) {
								$badC ++ if $valz[$m] =~ /^(5|9)$/;
							}
							else {
								$goodC ++ if $valz[$m] =~ /^(5|9)$/;
							}
						}
						if ($goodC < 5 and $badC >= 9) {
#							print "\t$LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
							print "\t$YW$name$N $LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n";
							$checkBad = 1; last;
						}
						else {
							#print "\t$LGN OKAY!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
						}
	###					print "\t\t$LGN YES$N\n" if $name eq "SEQ_77011";
					#}
				}
			}
#			print "\t$name checkbad = $checkBad\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746";
			
			if ($checkBad == 1) {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n";
###				print "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n" if $name eq "SEQ_77011";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
			elsif ($end > 100 + $edge1 and $beg < $edge2-100) {
###				print "something wrong\n" if $name eq "SEQ_77011";
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Used: $LGN$peak$N\n";
				push(@peak, "$beg-$end");
				push(@{$peak2{peak}}, $peak);
			}
			else {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Not : $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
		}
	}
	my $totalpeak = scalar(@peak);
	my @val2;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		$val2[$i] = $val;
		if ($val =~ /^(5|9)$/ and defined $nopeak{$i}) { # Peak Converted CpG or CH
			$val2[$i] = 4 if $val eq 5;
			$val2[$i] = 1 if $val eq 9;
		}
	}
	#die $print if $totalpeak > 1;
	$print .= "$name\t$totalpeak\n";
	#die "$print\n" if $totalpeak == 1;# or $name eq "SEQ_100022";
	return ($name, \@val2, $totalpeak, $peak2{peak});
}
__END__
# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

   # For nucleotide
# 10 = Nucleotide A
# 11 = Nucleotide C
# 12 = Nucleotide T
# 13 = Nucleotide G

