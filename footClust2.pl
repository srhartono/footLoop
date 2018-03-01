#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_g);
getopts("vd:n:g:");

#########
# BEGIN #
#########

BEGIN {
   my ($bedtools) = `bedtools --version`;
   my ($bowtie2) = `bowtie2 --version`;
   my ($bismark) = `bismark --version`;
   my ($bismark_genome_preparation) = `bismark_genome_preparation --version`;
	print "\n\n\e[1;33m ------------- BEGIN ------------ \e[0m\n";
   if (not defined $bedtools or $bedtools =~ /command not found/ or $bedtools =~ /bedtools v?([01].\d+|2\.0[0-9]|2\.1[0-6])/) {
      print "Please install bedtools at least version 2.17 before proceeding!\n";
      $bedtools = 0;
   }
   print "\n- \e[1;32m bedtools v2.17+ exists:\e[0m " . `which bedtools` if $bedtools ne 0;
   die if $bedtools eq 0;
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
#	print "\n\n\e[1;33m ------------ BEGIN ------------> \e[0m\n";
}
use myFootLib; use FAlite;

################
# ARGV Parsing #
###############


my $date = getDate();

my ($footPeakFolder) = ($opt_n);

# sanity check -n footPeakFolder
die "\nUsage: $YW$0$N $CY-n <footPeak's output folder (footPeak's -o)>$N\n\n" unless defined $opt_n and -d $opt_n;
($footPeakFolder) = getFullpath($footPeakFolder);
my $outDir = "$footPeakFolder/FOOTCLUST/";
makedir("$outDir/.CALL/");
die "Failed to create output directory $LCY$outDir/.CALL/$N!\n" unless -d "$outDir/.CALL";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footClust2.txt") or die "Failed to create outLog file $outDir/logFile_footClust2.txt: $!\n";
LOG($outLog, "\n\n$YW -------- PARSING LOG FILE -------- $N\n\n");

# get .fa file from footPeakFolder and copy
my ($faFile) = <$outDir/*.fa>;
LOG($outLog, date() . "\n$LRD Error$N: Cannot find any fasta file in $LCY$footPeakFolder$N\n\n") if not defined $faFile or not -e $faFile;
$faFile =~ s/\/+/\//g;
my %fa;
open (my $inFA, "<", $faFile) or die;
my $fasta = new FAlite($inFA);
while (my $entry = $fasta->nextEntry()) {
	my $def = $entry->def; $def =~ s/^>//; $def = uc($def);
	my $seq = $entry->seq;
	$fa{$def} = $seq;
}
close $inFA;

# get .local.bed peak files used for clustering
my $folder;
my %data;
open (my $inLog, "<", "$outDir/.0_LOG_FILESRAN") or DIELOG($outLog, date() . __LINE__ . "Failed to read from $outDir/.0_LOG_FILESRAN: $!\n");
while (my $line = <$inLog>) {
	chomp($line);
	if ($line =~ /^#FOLDER/) {
		($folder) = $line =~ /^#FOLDER=(.+)$/;
		next;
	}
	my ($pcb, $gene, $strand, $window, $thres, $type, $skip, $total_peak_all, $total_read_unique, $total_peak_used, $peaks_local_file, $peaks_file, $cluster_file) = split("\t", $line);
	$gene = uc($gene);
	$data{skip}{$gene} ++ if $skip =~ /Skip/;
	next if $skip =~ /Skip/;
	($cluster_file) =~ s/^CLUSTER_FILE=(.+)$/$1/; 
	DIELOG($outLog, date() . __LINE__ . "Failed to parse cluster file from $cluster_file\n") if $cluster_file =~ /\=/;
	($peaks_file) =~ s/^PEAK_FILE=(.+)$/$1/; 
	DIELOG($outLog, date() . __LINE__ . "Failed to parse peaks file from $peaks_file\n") if $peaks_file =~ /\=/;
	($peaks_local_file) =~ s/^PEAKS_LOCAL=(.+)$/$1/; 
	DIELOG($outLog, date() . __LINE__ . "Failed to parse peaks_local file from $peaks_local_file\n") if $peaks_local_file =~ /\=/;
	$strand = $cluster_file =~ /_Pos_/ ? "Pos" : $cluster_file =~ /_Neg_/ ? "Neg" : "Unk";
	my $type2 = (($strand eq "Pos" and $type =~ /^C/) or ($strand eq "Neg" and $type =~ /^G/)) ? "${LGN}GOOD$N" : "${LRD}WEIRD$N";
	LOG($outLog, date() . "$type2, pcb=$pcb, gene=$gene, strand=$strand, window=$window, thres=$thres, type=$type, skip=$skip, total_peak_all=$total_peak_all, total_read_unique=$total_read_unique, total_peak_used=$total_peak_used, peaks_local_file=$peaks_local_file, peaks_file=$peaks_file, cluster_file=$cluster_file\n");
	$data{good}{$gene} ++;
	%{$data{file}{$gene}{$cluster_file}} = (
		folder => $folder,
		pcb => $pcb,
		gene => $gene,
		strand => $strand,
		window => $window,
		thres => $thres,
		type => $type,
		total_peak_all => $total_peak_all,
		total_read_unique => $total_read_unique,
		total_peak_used => $total_peak_used,
		peaks_file => $peaks_file,
		peaks_local_file => $peaks_local_file
	);
}
my $skipped = 0;
foreach my $gene (sort keys %{$data{skip}}) {
	next if defined $data{good}{$gene};
	LOG($outLog, date() . "${LRD}Skipped $LCY$gene$N\n");
	$skipped ++;
}
LOG($outLog, "\n" . date() . "\tSkipped $LGN$skipped$N genes!\n\n");

my $cluster_count = -1;
foreach my $gene (sort keys %{$data{file}}) {
	foreach my $cluster_file (sort keys %{$data{file}{$gene}}) {
		$cluster_count ++;
		my $strand2 = $cluster_file =~ /_Pos_/ ? "Pos" : $cluster_file =~ /_Neg_/ ? "Neg" : "Unk";
		my ($thres, $window, $type) = $cluster_file =~ /^.*$gene\_$strand2\_(\d+\.?\d*)_(\d+\.?\d*)_(CG|CH|GH|GC)/;
#		die "thres = $thres,winow=$window, file=$cluster_file\n";
#		next if $cluster_file =~ /(CG|GC)/;
		my $folder = $data{file}{$gene}{$cluster_file}{folder};
		my $fasta = $folder . "/" . $cluster_file . ".fa";
		my $bedFile   = $folder . "/" . $cluster_file . ".bed";
		die "Fasta does not exist! ($fasta)\n" if not -e $fasta or -s $fasta == 0;
		die "Bed does not exist! ($bedFile)\n" if not -e $bedFile or -s $bedFile == 0;
		print "parsing bed of gene $gene bedFile $bedFile, fasta=$fasta\n";
		my ($BED) = parse_bed($bedFile, $fa{$gene});
		$BED = parse_fasta($BED, $fasta);
		$BED = shuffle_orig($BED);
		my %bed = %{$BED};
		my $want = "coor|gene|beg|end|cluster|total_peak|strand|pos|len|cpg|gc|skew";
		open (my $out2, ">", "$footPeakFolder/$cluster_file.kmer") or die "Cannot write to $LCY$footPeakFolder/$cluster_file$N: $!\n";
		foreach my $coor (sort {$bed{$a}{"1d_cluster"} cmp $bed{$b}{"1d_cluster"}} keys %bed) {
			print $out2 "coor\ttype" if $cluster_count == 0;
			foreach my $key (sort keys %{$bed{$coor}}) {
				next if $key !~ /($want)/;
				my ($header) = $key =~ /^\d+[a-z]_(\w+)$/;
					($header) = $key =~ /^.+_(\w+)$/ if not defined $header;
				if ($bed{$coor}{$key} =~ /^HASH/ and defined $bed{$coor}{$key}{shuf} and $bed{$coor}{$key}{shuf} ne "NA") {
					print $out2 "\t$header.orig\t$header.shuf\t$header.odds\t$header.pval" if $cluster_count == 0;
				}
				else {
					print $out2 "\t$header" if $cluster_count == 0;
				}
			}
			last;
		}
		print $out2 "\n";
		foreach my $coor (sort {$bed{$a}{"1d_cluster"} cmp $bed{$b}{"1d_cluster"}} keys %bed) {
			print $out2 "$coor\t$type";
			foreach my $key (sort keys %{$bed{$coor}}) {
				next if $key !~ /($want)/;
				if ($bed{$coor}{$key} =~ /^HASH/) {# and defined $bed{$coor}{$key}{shuf} and $bed{$coor}{$key}{shuf} ne "NA") {
					my $orig = $bed{$coor}{$key}{orig};
					my $shuf = defined ($bed{$coor}{$key}{shuf}) ? $bed{$coor}{$key}{shuf} : "NA";
					my $pval = defined ($bed{$coor}{$key}{pval}) ? $bed{$coor}{$key}{pval} : "NA";
					my $odds = defined ($bed{$coor}{$key}{odds}) ? $bed{$coor}{$key}{odds} : "NA";
					if ($bed{$coor}{$key}{shuf} ne "NA") {
						print $out2 "\t$orig\t$shuf\t$odds\t$pval";
					}
					else {
						print $out2 "\t$orig";
					}
				}
				else {
#				if ($key !~ /^2\w_/ or not defined $bed{$coor}{$key}{orig} or not defined $bed{$coor}{$key}{shuf} ) {
					print $out2 "\t$bed{$coor}{$key}";
				}
			}
			print $out2 "\n";
		}
		print "Done\n";
	}
}


sub parse_bed {
	my ($bedFile, $ref) = @_;
	my %bed;
	open (my $in, "<", $bedFile) or die;
	while (my $line = <$in>) {
		chomp($line);
		next if ($line =~ /^#/);
		my ($gene, $beg, $end, $cluster, $total_peak, $strand) = split("\t", $line);
		next if $end - $beg < 10;
		my $coor = "$gene:$beg-$end($strand)";
		my ($CLUST, $POS) = $cluster =~ /^(\d+)\.(.+)$/;
		if ($strand eq "-") {
			$CLUST = "END" if $cluster =~ /BEG/;
			$CLUST = "BEG" if $cluster =~ /END/;
		}
		%{$bed{$coor}} = (
			'1a_gene' => $gene,
			'1b_beg' => $beg,
			'1c_end' => $end,
			'1d_cluster' => $CLUST,
			'1e_total_peak' => $total_peak,
			'1f_strand' => $strand,
			'1g_pos' => $POS
		);
		my ($ref1) = $ref =~ /^(.{$beg})/;
		my $lenref = length($ref) - $end;
		my ($ref2) = $ref =~ /(.{$lenref})$/;
		print "gene $gene beg=$beg, end=$end, length = " . length($ref) . ", lenref= $lenref\n";
		$ref1 = "" if not defined $ref1;
		$ref2 = "" if $lenref < 1 or not defined $ref2;
		$ref1 .= "N$ref2";
		if ($strand eq "-" or $strand eq "Neg") {
			$ref1 = revcomp($ref1);
		}
		$bed{$coor}{'1h_ref'} = $ref1;
		# shuffle the ref1
#		my $times = 20;
#		my $shuf = shuffle_fasta($ref1, $end - $beg, $times);
		#$bed{$coor} = calculate_gcprofile($bed{$coor}, $ref1, 2, "shuf");
	}
	close $in;
	return (\%bed);
}

sub shuffle_orig {
	my ($BED) = @_;
	foreach my $coor (sort keys %{$BED}) {
		my $ref = $BED->{$coor}{'1h_ref'};
		my $lenseq = $BED->{$coor}{'2a_len'}{orig};
		die "coor = $coor\n" if not defined $lenseq;
		die "Undefined ref of gene $coor and lenref\n" if not defined $ref;

		$BED->{$coor} = shuffle_fasta($BED->{$coor}, $ref, $lenseq, 1000);
	}
	return $BED;
}

sub shuffle_fasta {
	my ($BED, $ref, $lenseq, $times) = @_;
#	$ref = substr($ref, 0, 250);
	my $lenref = length($ref);
	my $gcprof;
	for (my $i = 0; $i < $times; $i++) {
		my $ref2 = "";
		# if length of remaining ref seq is less than length of peak, then randomly take 100bp chunks until same length
		if ($lenref < $lenseq) {
			my $count = 0;
			my $add = int(0.4*$lenref);
			while (length($ref2) < $lenseq + $count) {
				my $randbeg = int(rand($lenref-$add));
				my $add2 = length($ref2) > $lenseq + $count - $add ? ($lenseq + $count - length($ref2)) : $add;
				$ref2 .= "N" . substr($ref, $randbeg, $add2);
			#	print "$count: add=$add, $randbeg + $add2, len = " . length($ref2) . " les than $lenseq + $count: $ref2\n";
				$count ++;
			}
			#print "len = " . length($ref2) . ", length seq = $lenseq + $count\n";
		}
		else {
			my $randbeg = int(rand($lenref-$lenseq));
			$ref2 = substr($ref, $randbeg, $lenseq);
			
		}
#		print "len = " . length($ref2) . ", length seq = $lenseq\n";
		$gcprof->[$i] = calculate_gcprofile($gcprof->[$i], $ref2, 2, "shuf");
	}
#		print "$i";
	foreach my $key (sort keys %{$gcprof->[0]}) {
		if ($key !~ /(cpg|gc|skew|kmer)/) {
			$BED->{$key}{shuf} = "NA";
			$BED->{$key}{pval} = "NA";
			$BED->{$key}{odds} = "NA";
			next;
		}
		my $orig = $BED->{$key}{orig}; die "key=$key,\n" if not defined $orig;
		my ($nge, $nle, $ntot) = (0,0,scalar(@{$gcprof}));
		my @temp;
		for (my $i = 0; $i < @{$gcprof}; $i++) {
			my $shuf = $gcprof->[$i]{$key}{shuf};
			$nge ++ if $orig >= $shuf;
			$nle ++ if $orig <= $shuf;
#			print "\t$gcprof->[$i]{$key}{shuf}" if $key =~ /(cpg|gc|skew)/;
			$temp[$i] = $gcprof->[$i]{$key}{shuf} if $key =~ /(cpg|gc|skew|kmer)/;
		}
		$BED->{$key}{pval} = $nge > $nle ? myformat(($nle+1)/($ntot+1)) : myformat(($nge+1)/($ntot+1));
		$BED->{$key}{shuf} = myformat(tmm(@temp));
		my $shuf = $BED->{$key}{shuf};
		my $pval = $BED->{$key}{pval};
		$BED->{$key}{odds} = ($orig > 0 and $shuf > 0) ? (0.1+$orig) / (0.1+$shuf) : ($orig < 0 and $shuf < 0) ? (abs($shuf)+0.1) / (abs($orig)+0.1) : $orig < 0 ? -1 * (abs($orig)+0.1) / 0.1 : (abs($orig+0.1)) / 0.1;
		$BED->{$key}{odds} = myformat($BED->{$key}{odds});
		my $odds = $BED->{$key}{odds};
#		print "$key $odds: $orig vs. $shuf, p=$pval\n";
	}
	return $BED;
}

sub parse_fasta {
	my ($bed, $fastaFile) = @_;
	open (my $in, "<", $fastaFile) or die;
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def; $def =~ s/^>//;
		my $seq = $entry->seq;
		next if length($seq) < 10;
		$bed->{$def} = calculate_gcprofile($bed->{$def}, $seq, 2, "orig");
#		print "def=$def lenseq=$bed->{$def}{'2a_len'}{orig}\n";
	}
	close $in;
	return $bed;
}

sub calculate_gcprofile {
	my ($bed, $seq, $number, $type) = @_;
	$seq = uc($seq);
	my ($G) = $seq =~ tr/G/G/;
	my ($C) = $seq =~ tr/C/C/;
	my ($A) = $seq =~ tr/A/A/;
	my ($T) = $seq =~ tr/T/T/;
	my ($N) = $seq =~ tr/N/N/;
	my $CG = 0;
	while ($seq =~ /CG/g) {
		$CG ++;
	}
	my $len = length($seq) - length($N);
	my $gc = int(100 * ($G+$C) / $len+0.5)/100;
	my $skew  = $gc == 0 ? 0 : int(100 * (abs($G-$C)+1)  / ($G+$C+1)  +0.5)/100;
		$skew *= -1 if $C > $G;
	my $skew2 = $gc == 0 ? 0 : int(100 * (abs($G-$C)+10) / ($G+$C+10) + 0.5)/100;
		$skew2 *= -1 if $C > $G;
	my $cpg = int(100 * ($CG+1) / (($C+1) * ($G+1)) * $len + 0.5)/100;
	$bed->{$number . 'a_len'}{$type} = $len;
	$bed->{$number . 'b_cpg'}{$type} = $cpg;
	$bed->{$number . 'c_gc'}{$type} = $gc;
	$bed->{$number . 'd_skew'}{$type} = $skew;
	$bed->{$number . 'e_skew2'}{$type} = $skew2;
	$bed->{$number . 'f_acgt'}{$type} = "CG=$CG, C=$C, G=$G, A=$A, T=$T, len=$len";
	return $bed;
}


sub calculate_kmer {
	my ($seq) = @_;
	

}

__END__
my @local_peak_files = <$footPeakFolder/PEAKS_LOCAL/*.local.bed>;
LOG($outLog, date() . "\nError: cannot find any .local.bed peak files in $LCY$footPeakFolder$N\n\n") and die if not -d "$footPeakFolder/PEAKS_LOCAL/" or @local_peak_files == 0;

# Log initialization info
my $init = "\n\n$YW ------------- INIT ------------- $N\n" . "\n
Date                       = $date;
Command                    = $0 -d $dist -n $footPeakFolder
${LGN}FootPeakFolder    $N = $footPeakFolder
${LGN}Fasta File        $N = $faFile
${LGN}Inbetween Min Dist$N = $dist\n\n
$YW -------- RUN (" . scalar(@local_peak_files) . " files) --------- $N\n\n";
LOG($outLog, $init);


####################
# Processing Input #
####################
my $files_log = "#FOLDER=$footPeakFolder\n";
my $input1_count = -1;
foreach my $input1 (sort @local_peak_files) {
	$input1_count ++;
#DEBUG
#	next if $input1 !~ /CALM3/;	

	# remove double // from folder
	$input1 =~ s/[\/]+/\//g;
	my ($folder1, $fileName1) = getFilename($input1, "folder");
	my ($fullName1) = getFilename($input1, "full");

	# get gene and strand from file name
	my ($gene, $strand, $window, $thres, $type) = $fileName1 =~ /^(.+)\_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(GH|GC|CG|CH).PEAK/;
	LOG($outLog, date() . "Cannot parse gene name from file=$LCY$input1$N\n") unless defined $gene and defined $strand and defined $window and defined $thres and defined $type;
	$thres *= 100 if $thres < 1;
	$gene   = uc($gene);
	$strand = $strand eq "Neg" ? "-" : "+";
	my ($pcb) = $footPeakFolder =~ /(PCB\d+)/; $pcb = "PCBUNK" if not defined $pcb;
	# check peak file. Skip if there's less than 10 peaks
	my ($total_line) = `wc -l $input1` =~ /^(\d+)/;
	if ($total_line < 10) {
		LOG($outLog, "\n--------------\n\n$LGN$input1_count$N. ${LRD}NEXTED $N: $input1 ${LCY}because total peaks is less than 10 $N($LGN$total_line$N)\n\n");
		$files_log .= "$pcb\t$gene\t$strand\t$window\t$thres\t$type\t${LRD}Skip$N\ttotal_peak_all=$total_line\ttotal_read_unique=-1\ttotal_peak_used=-1\tPEAKS_LOCAL=PEAKS_LOCAL/$fullName1\tPEAK_FILE=NA\n";
		next;
	}
	$files_log .= "$pcb\t$gene\t$strand\t$window\t$thres\t$type\t${LGN}Ran$N\ttotal_peak_all=$total_line";

	# Actual parse of peak file
	my ($linecount, $total_peak_used, $total_peak_all, %data) = (0,0,0);
	open (my $in1, "sort -k1,1 -k2,2n -k3,3n $input1|") or LOG($outLog, date() . "Cannot read from $input1: $!\n") and die;
	LOG($outLog, "\n--------------\n\n$LGN$input1_count$N. ${LCY}RUNNING$N: $input1\n");
	LOG($outLog, date() . "$LCY\tInfo$N: pcb=$pcb,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type\n");
	LOG($outLog, date() . "$LCY\tExample Lines$N:\n");
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my ($read, $beg, $end) = split("\t", $line);
		my ($num) = $read =~ /^.+\/(\d+)\/ccs/; 
		LOG($outLog, date() . "$LRD\tERROR$N:$LCY Read must end in this format$N: <anything>/<hole number>/ccs\n\n$read\n\n") and die if not defined $num;
		my $check = 0;
		$total_peak_all ++;

		# reads with multiple peaks separated by less than distance (250) is merged together
		if (defined $data{$num}) {
			for (my $i = 0; $i < @{$data{$num}}; $i++) {
				my $beg2 = $data{$num}[$i][0];
				my $end2 = $data{$num}[$i][1];
				if ($beg < $end2 + $dist) {
					$data{$num}[$i][1] = $end;
					$check = 1;
					LOG($outLog, date() . "$LGN\tline=$linecount$N: readnum=$num MERGED beg2=$beg2 end2=$end2 with beg=$beg end=$end into beg3=$beg2 end3=$end\n") if $linecount < 5;
					last;
				}
			}
		}

		# if 1 peak or peak is far then create new peak
		if ($check == 0) {
			push(@{$data{$num}}, [$beg,$end,$read,$gene]);
			$total_peak_used ++;
			LOG($outLog, date() . "$LGN\tline=$linecount$N: readnum=$num beg=$beg end=$end\n") if $linecount < 5;
		}
	}
	close $in1;

	# Put info into files log for next script	
	my $total_read_unique = (keys %data);
	my ($fileNameShort) = $fullName1 =~ /^(.+_(GH|CH|GC|CG))/;
	die "Cannot get filenameshort of $input1 ($LCY$fullName1$N)\n" unless defined $fileNameShort;
	my $peakFile = "$footPeakFolder/.CALL/$fileNameShort.PEAK.out";
	die "Cannot find PEAK file of $input1 ($LCY$peakFile$N)\n" unless -e $peakFile;
	$files_log .= "\ttotal_read_unique=$total_read_unique\ttotal_peak_used=$total_peak_used\tPEAKS_LOCAL=PEAKS_LOCAL/$fullName1\tPEAK_FILE=$peakFile\tCLUSTER_FILE=FOOTCLUST/.TEMP/$fullName1.clust\n";

	# Calculate average beg/mid/end and write a temp bed file for each peak
	LOG($outLog, date() . "$LCY\tInfo 2$N: total_peak_all=$LGN$total_peak_all$N,total_read_unique=$LGN$total_read_unique$N,total_peak_used=$LGN$total_peak_used$N\n");
	open (my $out1, ">", "$outDir/.TEMP/$fullName1.temp") or DIELOG($outLog, date() . __LINE__ . "\tCannot write to $LCY$fileName1.temp$N: $!\n");
	print $out1 "id\tbeg\tend\n";
	foreach my $num (sort keys %data) {
		for (my $i = 0; $i < @{$data{$num}}; $i++) {
			my $beg = $data{$num}[$i][0];
			my $end = $data{$num}[$i][1];
			my $read = $data{$num}[$i][2];
			print $out1 "$num.$i\t$beg\t$end\n";
		}
	}
	close $out1;

	# Write and run R script
	LOG($outLog, date() . "$LCY\tRunning R script$N $outDir/.TEMP/$fullName1.temp.R\n");
	my $Rscript = "
	
	set.seed(420)
	setwd(\"$outDir\");
	library(ggplot2)
	df = read.table(\"$outDir/.TEMP/$fullName1.temp\",row.names=1,header=T,sep=\"\\t\")
	dm = kmeans(df,5,nstart=20)
	df\$cluster = dm\$cluster
	df = df[order(df\$cluster, df[,1], df[,2]),]
	df\$y = seq(1,dim(df)[1])
	df\$ymax = seq(2,dim(df)[1]+1)
	colnames(df) = c(\"x\",\"xmax\",\"clust\",\"y\",\"ymax\")
	df2 = as.data.frame(aggregate(df[,c(1,2)],by=list(df\$clust),function(x)mean(x,trim=0.05)))
	colnames(df2) = c(\"clust\",\"x2\",\"y2\")
	df2 = df2[order(df2\$x2, df2\$y2),]
	df2\$clust2 = seq(1,dim(df2)[1])
	df\$id = rownames(df)
	df = as.data.frame(merge(df,df2[,c(1,4)]))
	df\$clust = df\$clust2
	df = df[order(df\$clust, df\$x, df\$xmax),]
	df\$y = seq(1,dim(df)[1])
	df\$ymax = seq(2,dim(df)[1]+1)
	df[,c(1,6)] = df[,c(6,1)]
	colnames(df)[c(1,6)] = colnames(df)[c(6,1)]
	df = df[,-7]
	png(\"$outDir/PNG/$fullName1.clust.png\",height=(dim(df)[1])*5,width=max(df\$xmax)+10)
	ggplot(df, aes(x,y)) + geom_rect(aes(xmin=x,ymin=y,xmax=xmax,ymax=ymax,fill=as.factor(clust))) + theme_bw() + 
	theme(panel.grid=element_blank()) + coord_fixed(ratio=10)
	dev.off()
	write.table(df,\"$outDir/.TEMP/$fullName1.clust\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
	";
	open (my $outR, ">", "$outDir/.TEMP/$fullName1.temp.R") or DIELOG($outLog, date() . __LINE__ . "Failed to write to $outDir/.TEMP/$fullName1.temp.R: $!\n");
	print $outR $Rscript;
	system("R --vanilla --no-save < $outDir/.TEMP/$fullName1.temp.R > $outDir/.TEMP/$fullName1.temp.R.log 2>&1");
	close $outR;
	my @Rlog = `tail -n 1 $outDir/.TEMP/$fullName1.temp.R.log`;
	LOG($outLog, date() . "$LCY\tR logs$N:\n\t\t\t\t" . join("\n\t\t\t\t", @Rlog) . "\n");
	
	# process clust
	LOG($outLog, date() . "$LCY\tCreating fasta file for each cluster$N $outDir/.TEMP/$fullName1.clust.fa\n");
	my %cl;
	open (my $in2, "<", "$outDir/.TEMP/$fullName1.clust") or DIELOG($outLog, date() . __LINE__ . "\tFailed to read from $outDir/.TEMP/$fullName1.clust: $!\n");
	while (my $line = <$in2>) {
		chomp($line);
		next if $line =~ /(clust|xmax|ymax)/;
		my ($num, $beg, $end, $y, $y2, $clust) = split("\t", $line);
		my ($ind) = "";
		($num, $ind) = $num =~ /^(\d+)\.(\d+)$/ if $num =~ /^\d+\.\d+$/;
		LOG($outLog, date() . "Cannot parse num from line=$line\n") if not defined $num;
		my ($mid) = int(($end + $beg)/2+0.5);
		push(@{$cl{$clust}{beg}}, $beg);
		push(@{$cl{$clust}{end}}, $end);
	}
	close $in2;
	
	# process clust: get average beg/end point
	open (my $out2, ">", "$outDir/.TEMP/$fullName1.clust.bed") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outDir/.TEMP/$fullName1.clust.bed: $!\n");
	print $out2 "#gene\tbeg\tend\tcluster\ttotal_peak\tstrand\n";
	foreach my $clust (sort {$a <=> $b} keys %cl) {
		my $beg = int(tmm(@{$cl{$clust}{beg}})+0.5);
		my $end = int(tmm(@{$cl{$clust}{end}})+0.5);
		my $total = @{$cl{$clust}{beg}};
		print $out2 "$gene\t$beg\t$end\t$clust\t$total\t$strand\n";
	}
	close $out2;
	system("fastaFromBed -fi $faFile -bed $outDir/.TEMP/$fullName1.clust.bed -fo $outDir/.TEMP/$fullName1.clust.fa -s");

	LOG($outLog, date() . "$LCY\tAnalyzing Sequence Content from$N $outDir/.TEMP/$fullName1.clust.fa\n");
	open (my $faIn, "$outDir/.TEMP/$fullName1.clust.fa") or DIELOG($outLog, "Failed to read from $outDir/.TEMP/$fullName1.clust.fa: $!\n");
	my $fasta = new FAlite($faIn); $linecount = 0;
	while (my $entry = $fasta->nextEntry()) {
		$linecount ++;
		my $def = uc($entry->def); $def =~ s/^>//;
		my $seq = uc($entry->seq());
		print ">$def\n$seq\n\n" if $linecount == 1;
	}
}	
open (my $outz, ">", "$outDir/.0_LOG_FILESRAN") or DIELOG($outLog, "Failed tow rite to $outDir/.0_LOG_FILESRAN:$!\n");
print $outz $files_log;
close $outz;
LOG($outLog, "\n$YW ----------- FILES RAN ---------- $N\n\n" . date() . "$LCY\tSummary of run$N: $outDir/.0_LOG_FILESRAN\n\n$files_log\n\n")
	
	
	
	
	
	
	
__END__
	close $out1;
	
	#df = read.table("CALM3_Pos_20_0.65_CH.PEAK.local.bed.temp",row.names=1,header=F,sep="\t")
	
	set.seed(420)
	library(ggplot2)
	df = read.table("CALM3_Pos_20_0.65_CH.PEAK.local.bed.temp",row.names=1,header=T,sep="\t")
	dm = kmeans(df,6,nstart=20)
	df$cluster = dm$cluster
	df = df[order(df$cluster, df[,1], df[,2]),]
	df$y = seq(1,dim(df)[1])
	df$ymax = seq(2,dim(df)[1]+1)
	colnames(df) = c("x","xmax","clust","y","ymax")
	df2 = as.data.frame(aggregate(df[,c(1,2)],by=list(df$clust),function(x)mean(x,trim=0.05)))
	colnames(df2) = c("clust","x2","y2")
	df2 = df2[order(df2$x2, df2$y2),]
	df2$clust2 = seq(1,dim(df2)[1])
	df$id = rownames(df)
	df = as.data.frame(merge(df,df2[,c(1,4)]))
	df$clust = df$clust2
	df = df[order(df$clust, df$x, df$xmax),]
	df$y = seq(1,dim(df)[1])
	df$ymax = seq(2,dim(df)[1]+1)
	df[,c(1,6)] = df[,c(6,1)]
	colnames(df)[c(1,6)] = colnames(df)[c(6,1)]
	df = df[,-7]
	png("CALM3_Pos_20_0.65_CH.PEAK.local.bed.clust.png",height=(dim(df)[1])*5,width=max(df$xmax)+10)
	ggplot(df, aes(x,y)) + geom_rect(aes(xmin=x,ymin=y,xmax=xmax,ymax=ymax,fill=as.factor(clust))) + theme_bw() + 
	coord_fixed(ratio=10)
	dev.off()
	
