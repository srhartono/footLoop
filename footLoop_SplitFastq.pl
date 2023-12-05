#!/usr/bin/perl
# version 231204

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_i $opt_n $opt_o $opt_l);
getopts("vi:n:o:l:");

my ($input1, $outDir, $lenthreshold) = ($opt_i, $opt_o, $opt_l);

my $usage = "\nUsage: $YW$0$N -i $LCY<fastq>$N -o $LCY<output DIRECTORY>$N -n $LGN<number of reads per file [1000000]>$N\n\n";
die $usage unless defined($opt_i) and defined($opt_o) and -e $opt_i;

# splitby
my $splitBy = defined $opt_n ? $opt_n : 1000000;
if ($splitBy !~ /^\d+$/ or $splitBy eq 0) {
	die date() . $usage . "\n$LRD!!! Fatal error: -n ${YW}have to be whole positive integer number!$N currently:$LGN$splitBy$N\n\n";
}

# Create Directory
if (not -d $outDir) {
	my $system_mkdir = system("mkdir -p $outDir");
	if ($system_mkdir != 0) {
		print date() . "$LRD!!!$N FATAL ERROR!!! Failed to create directory $CY$outDir$N $YW(make sure the previous directory exists!)$N\n\n";
		die "\n";
	}
	else {
		print date() . "\nCreated directory $CY$outDir$N\n";
	}
}


my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my ($readNumber, $printExample, $meanSeqLength, $totalRead) = (0,0);

my $in1;
if ($input1 =~ /.gz$/) {
	open ($in1, "zcat $input1|") or die date() . "Cannot read from $input1: $!\n";
}
else {
	open ($in1, "<", $input1) or die date() . "Cannot read from $input1: $!\n";
}

my $outFile = "$outDir/$fileName1"; #.0.part.gz";
my $out1;
open ($out1, "| gzip > $outFile.0.part.gz") or die date() . "Cannot write to $outFile.0.part.gz: $!\n";
print date() . "${YW}First sequence found:$N\n";

my $linecount = 0;
while (my $line = <$in1>) {
	chomp($line);
	$linecount ++;

	if ($readNumber % $splitBy == 0) {
		my $number = int($readNumber / $splitBy);
		print date() . "$input1: Done $readNumber\n";
		open ($out1, "| gzip > $outFile.$number.part.gz") or die date() . "Cannot write to $outFile.$number.part.gz: $!\n";
	}

	# name
	print $out1 "$line\n";
	print date() . "${LPR}Name$N\t$line\n" if $printExample == 0;
	my $name = $line; $name =~ s/^\@//;

	$line = <$in1>; chomp($line); $linecount ++;

	# sequence
	print $out1 "$line\n";
	my $seq      = $line;
	my $seqprint = $seq;
		$seqprint = substr($seq, 0, 30) . "... (" . length($seq) . " bp)" if length($seq) > 30;
	$meanSeqLength += length($seq); $totalRead ++;
	print date() . "${LPR}Seq$N\t$seqprint\n" if $printExample == 0;

	$line = <$in1>; chomp($line); $linecount ++;

	# plus
	print $out1 "$line\n";
	print date() . "${LPR}Name2$N\t$line\n" if $printExample == 0;
	my $plus = $line; $plus =~ s/^\+//; my $name2 = $line;

	$line = <$in1>; chomp($line); $linecount ++;

	#qual
	print $out1 "$line\n";
	my $qual = $line;
	my $qualprint = substr($qual, 0, 30) . "... (" . length($qual) . " bp)" if length($qual) > 30;
	print date() . "${LPR}Qual$N\t$qualprint\n" if $printExample == 0;

	if ($name2 !~ /^\+/ and $name !~ /$plus/i) {
		print "\n" . date() . "$input1: linecount $linecount: DIED AT READ $readNumber because line 3 of each read HAS to be either a${CY}+$N or its own name, which in this case it isn't:\n";
		print "  1: NAME:$name\n";
		print "  2: SEQ :$seq\n";
		print "  3: PLUS:$name2\n";
		print "  4: QUAL:$qual\n" . "\n";
		die "Corrupted fastq file due to above!\n\n";
	}

	$printExample = 1;
	$readNumber ++;
}
close $in1;
close $out1;

$meanSeqLength = $totalRead == 0 ? 0 : int($meanSeqLength / $totalRead * 100 + 0.5)/100;
my $totalFiles = $totalRead % $splitBy == 0 ? int($totalRead / $splitBy) : int($totalRead / $splitBy) + 1;
print date() . "$YW$input1$N: Total_read=$LGN$totalRead$N, Split_by=$LGN$splitBy$N, Total_files=$LGN$totalFiles$N, LastReadNumber=$LGN$readNumber$N, Mean_read_length=$LGN$meanSeqLength$N\n\n";
print "$totalRead,$splitBy,$totalFiles,$meanSeqLength$N\n";
