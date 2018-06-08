#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G);
getopts("vd:n:G:");

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
my ($genewant) = $opt_G if defined $opt_G;

# sanity check -n footPeakFolder
die "\nUsage: $YW$0$LGN [Optional: -G (genewant)]$N $CY-n <footPeak's output folder (footPeak's -o)>$N\n\n" unless defined $opt_n and -d $opt_n;
($footPeakFolder) = getFullpath($footPeakFolder);
my $outDir = "$footPeakFolder/FOOTCLUST/";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footPeak_kmer.TEMP.txt") or die "Failed to create outLog file $outDir/logFile_footPeak_kmer.txt: $!\n";
LOG($outLog, " \n\n$YW -------- PARSING LOG FILE -------- $N\n\n");

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

my $label = "";
if (-e "$footPeakFolder/.LABEL") {
	($label) = `cat $footPeakFolder/.LABEL`;
	chomp($label);
}
else {
	DIELOG($outLog, "Failed to parse label from .LABEL in $footPeakFolder/.LABEL\n");
}

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
	my ($label2, $gene, $strand, $window, $thres, $type, $skip, $total_peak_all, $total_read_unique, $total_peak_used, $peaks_local_file, $peaks_file, $cluster_file) = split("\t", $line);
   if (defined $opt_G and $line !~ /$opt_G/i) {
      LOG($outLog, date() . " Skipped $LCY$gene$N as it doesn't contain $LGN-G $opt_G$N\n","NA");
      next;
   }
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
	my $type2 = (($strand eq "Pos" and $type =~ /^C/) or ($strand eq "Neg" and $type =~ /^G/)) ? "${LGN}CORRECT_CONVERSION$N" : "${YW}REVERSED_CONVERSION$N";
	LOG($outLog, date() . "$type2, label=$label, gene=$gene, strand=$strand, window=$window, thres=$thres, type=$type, skip=$skip, total_peak_all=$total_peak_all, total_read_unique=$total_read_unique, total_peak_used=$total_peak_used, peaks_local_file=$peaks_local_file, peaks_file=$peaks_file, cluster_file=$cluster_file\n");
	$data{good}{$gene} ++;
	%{$data{file}{$gene}{$cluster_file}} = (
		folder => $folder,
		label => $label,
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
		my %print;
		$cluster_count ++;
		LOG($outLog, date() . " $YW$cluster_count$N. Parsing $LCY$cluster_file$N\n");
		my ($cluster_file_name) = getFilename($cluster_file, "full");
		my ($label2, $gene2, $strand, $window, $thres, $type, $barcode, $plasmid, $desc) = parseName($cluster_file_name);
		LOG($outLog, "Using label=$label2. Inconsistent label in filename $LCY$cluster_file_name$N\nLabel from $footPeakFolder/.LABEL: $label\nBut from fileName: $label2\n\n") if $label2 ne $label;
		$label = $label2;

		my $folder = $data{file}{$gene}{$cluster_file}{folder};
		my $fasta = $folder . "/" . $cluster_file . ".fa";
		my $bedFile   = $folder . "/" . $cluster_file . ".bed";
		if (not -e $fasta or -s $fasta == 0) {LOG($outLog, date() . " Fasta does not exist! ($fasta)\n"); next;}
		if (not -e $bedFile or -s $bedFile == 0) {LOG($outLog, date() . " Bed does not exist! ($bedFile)\n"); next;}
		LOG($outLog, date() . "\t- gene=$LGN$gene$N\n" . date() . "\t- bedFile=$LCY$bedFile$N\n" . date() . "\t- fasta=$LCY$fasta$N\n");
		my ($BED) = parse_bed($bedFile, $fa{$gene}, $outLog); 
		if (not defined $BED) {LOG($outLog, date() . " BED cannot be parsed from $bedFile so skipped!\n"); next;}
		$BED = parse_fasta($BED, $fasta, $outLog); next if not defined $BED;
		$BED = shuffle_orig($BED, $outLog);
		my %bed = %{$BED};
		my $want = "coor|gene|beg|end|cluster|total_peak|strand|pos|len|cpg|gc|skew2|ggg";
		my $kmerFile = $cluster_file; $kmerFile =~ s/.TEMP\/?//; 
		$kmerFile = "$footPeakFolder/$kmerFile.kmer";
		my $headerPrint = "";
		open (my $out2, ">", "$kmerFile") or print "Cannot write to $LCY$footPeakFolder/$cluster_file$N: $!\n" and next;
		my $keyz;
		foreach my $coor (keys %bed) {
			foreach my $key (keys %{$bed{$coor}}) {
				$keyz = "1d_cluster" if defined $bed{$coor}{"1d_cluster"};
				$keyz = $key if $key =~ /_cluster$/ and not defined $bed{$coor}{"1d_cluster"};
				last if defined $keyz;
			}
			last if defined $keyz;
		}
		if (not defined $keyz) {
			foreach my $coor (keys %bed) {
				foreach my $key (keys %{$bed{$coor}}) {
					$keyz = $key if $key =~ /_cluster$/;
					last if defined $keyz;
				}
				last if defined $keyz;
			}
		}
		foreach my $coor (sort {$bed{$a}{$keyz} cmp $bed{$b}{$keyz}} keys %bed) {
			my $typez = $coor =~ /BEG/ ? 1 : $coor =~ /MID/ ? 2 : $coor =~ /END/ ? 3 : $coor =~ /WHOLE/ ? 4 : 5;
			$headerPrint .= "coor\ttype";# if $cluster_count == 0;
			foreach my $key (sort keys %{$bed{$coor}}) {
				next if $key !~ /($want)/;
				my ($header) = $key =~ /^\d+[a-zA-Z]?_(\w+)$/;
					($header) = $key =~ /^.+_(\w+)$/ if not defined $header;
				if ($bed{$coor}{$key} =~ /^HASH/ and defined $bed{$coor}{$key}{shuf} and $bed{$coor}{$key}{shuf} ne "NA") {
					$headerPrint .= "\t$header.orig\t$header.shuf\t$header.odds\t$header.pval";
				}
				else {
					$headerPrint .= "\t$header";
				}
			}
			$headerPrint .= "\n";
			last;
		}
		foreach my $coor (sort {$bed{$a}{$keyz} cmp $bed{$b}{$keyz}} keys %bed) {
			my $typez = $coor =~ /BEG/ ? 1 : $coor =~ /MID/ ? 2 : $coor =~ /END/ ? 3 : $coor =~ /WHOLE/ ? 4 : 5;
			$print{$typez}{$coor} .= "$coor\t$type";
#			print "\tcoor=$coor typez=$typez\n";
			foreach my $key (sort keys %{$bed{$coor}}) {
#				print "\t\tkey=$key want=$want $LRD nexted$N\n" if $key !~ /($want)/;
				next if $key !~ /($want)/;
				if ($bed{$coor}{$key} =~ /^HASH/) {# and defined $bed{$coor}{$key}{shuf} and $bed{$coor}{$key}{shuf} ne "NA") {
					my $orig = $bed{$coor}{$key}{orig};
					my $shuf = defined ($bed{$coor}{$key}{shuf}) ? $bed{$coor}{$key}{shuf} : "NA";
					my $pval = defined ($bed{$coor}{$key}{pval}) ? $bed{$coor}{$key}{pval} : "NA";
					my $odds = defined ($bed{$coor}{$key}{odds}) ? $bed{$coor}{$key}{odds} : "NA";
					$orig = $orig eq "NA" ? "NA" : ($orig > 0.01 || $orig < -0.01) ? int($orig * 1000)/1000 : ($orig > 0.0001 || $orig < -0.0001) ? int($orig * 100000)/100000 : $orig;
					$shuf = $shuf eq "NA" ? "NA" : ($shuf > 0.01 || $shuf < -0.01) ? int($shuf * 1000)/1000 : ($shuf > 0.0001 || $shuf < -0.0001) ? int($shuf * 100000)/100000 : $shuf;
					$pval = $pval eq "NA" ? "NA" : ($pval > 0.01 || $pval < -0.01) ? int($pval * 1000)/1000 : ($pval > 0.0001 || $pval < -0.0001) ? int($pval * 100000)/100000 : $pval;
					$odds = $odds eq "NA" ? "NA" : ($odds > 0.01 || $odds < -0.01) ? int($odds * 1000)/1000 : ($odds > 0.0001 || $odds < -0.0001) ? int($odds * 100000)/100000 : $odds;
					if ($bed{$coor}{$key}{shuf} ne "NA") {
						$print{$typez}{$coor} .= "\t$orig\t$shuf\t$odds\t$pval";
#						print "\t\tkey=$key 1 orig=$orig shuf=$shuf\n";
					}
					else {
						$print{$typez}{$coor} .= "\t$orig";
#						print "\t\tkey=$key 2 orig=$orig shuf=NA\n";
					}
				}
				else {
#				if ($key !~ /^2\w_/ or not defined $bed{$coor}{$key}{orig} or not defined $bed{$coor}{$key}{shuf} ) {
#					print "\t\tkey=$key value=$bed{$coor}{$key}\n";
					$print{$typez}{$coor} .= "\t$bed{$coor}{$key}";
				}
			}
			$print{$typez}{$coor} .= "\n";
		}
		print $out2 "\#$headerPrint";
		foreach my $typez (sort {$a <=> $b} keys %print) {
			foreach my $coor (sort keys %{$print{$typez}}) {
				print $out2 $print{$typez}{$coor};
			}
		}
		LOG($outLog, "Done $cluster_file\n");
	}
}


sub parse_bed {
	my ($bedFile, $ref, $outLog) = @_;
	my %bed;
	my $linecount = 0;
	open (my $in, "<", $bedFile) or print "footPeak_kmer.pl: Cannot read from $bedFile: $!\n" and return;
	while (my $line = <$in>) {
		chomp($line);
		next if ($line =~ /^#/);
		my ($gene, $beg, $end, $cluster, $total_peak, $strand) = split("\t", $line);
		if (defined $opt_G and $gene !~ /$opt_G/i) {
			next;
		}
		next if $end - $beg < 10;
		$linecount ++;
		my $coor = "$gene:$beg-$end($strand)";
		my ($CLUST, $POS) = $cluster =~ /^(\d+)\.(.+)$/; 
		# $CLUST is cluster number of the cluster
		# $POS is the position type of cluster (BEG/MID/END/WHOLE)
		# $cluster contain the cluster and position type info (e.g. 1.WHOLE or 2.BEG)
		# e.g. $cluster is 17.BEG, then $CLUST = 17 and $POS = BEG
		if ($strand eq "-") {
			# BUG:
			$POS = "END" if $cluster =~ /BEG/;
			$POS = "BEG" if $cluster =~ /END/;
			#$CLUST = "END" if $cluster =~ /BEG/; #<- BUG
			#$CLUST = "BEG" if $cluster =~ /END/; #<- BUG
			# I should've put $POS instead of $CLUST, as $CLUST hold cluster number.
			# Otherwise, $CLUST will be changed from cluster number (e.g. 17) to type of cluster (e.g. END)
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
#		print "clster=$cluster beg=$beg end=$end total=$total_peak strand=$strand\n";
		foreach my $coor (sort keys %bed) {
			foreach my $key (sort keys %{$bed{$coor}}) {
#				print "\t\tkey=$key value=$bed{$coor}{$key}\n" if $key ne "1h_ref";
			}
		}
		my ($ref1) = $ref =~ /^(.{$beg})/;
		my $lenref = length($ref) - $end;
		my ($ref2) = $ref =~ /(.{$lenref})$/;
		LOG($outLog, date() . " \tbedfile line #$LGN$linecount$N: gene $LCY$gene$N beg=$LCY$beg$N, end=$LCY$end$N, length = " . length($ref) . ", lenref= $LCY$lenref$N\n","NA") if $linecount < 10;
		$ref1 = "" if not defined $ref1;
		$ref2 = "" if $lenref < 1 or not defined $ref2;
		$ref1 .= "N$ref2";
		if ($strand eq "-" or $strand eq "Neg") {
			$ref1 = revcomp($ref1);
		}
		$bed{$coor}{'1h_ref'} = $ref1;
	}
	close $in;
	return (\%bed);
}

sub shuffle_orig {
	my ($BED, $outLog) = @_;
	foreach my $coor (sort keys %{$BED}) {
		my $ref = $BED->{$coor}{'1h_ref'};
		my $lenseq = $BED->{$coor}{'2a_len'}{orig};
		DIELOG($outLog, date() . " Lenseq isn't defined coor = $coor\n") if not defined $lenseq;
		DIELOG($outLog, date() . " Lenseq isn't defined coor = $coor\n") if not defined $lenseq;
		DIELOG($outLog, date() . " Undefined ref of gene $coor and lenref\n") if not defined $ref;

		$BED->{$coor} = shuffle_fasta($BED->{$coor}, $ref, $lenseq, 1000);
	}
	return $BED;
}

sub shuffle_fasta {
	my ($BED, $ref, $lenseq, $times) = @_;
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
				$count ++;
			}
		}
		else {
			my $randbeg = int(rand($lenref-$lenseq));
			$ref2 = substr($ref, $randbeg, $lenseq);
			
		}
		$gcprof->[$i] = calculate_gcprofile($gcprof->[$i], $ref2, 2, "shuf");
	}
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
			$temp[$i] = $gcprof->[$i]{$key}{shuf} if $key =~ /(cpg|gc|skew|kmer)/;
		}
		$BED->{$key}{pval} = $nge > $nle ? myformat(($nle+1)/($ntot+1)) : myformat(($nge+1)/($ntot+1));
		$BED->{$key}{shuf} = myformat(tmm(@temp));
		my $shuf = $BED->{$key}{shuf};
		my $pval = $BED->{$key}{pval};
		$BED->{$key}{odds} = ($orig > 0 and $shuf > 0) ? (0.1+$orig) / (0.1+$shuf) : ($orig < 0 and $shuf < 0) ? (abs($shuf)+0.1) / (abs($orig)+0.1) : $orig < 0 ? -1 * (abs($orig)+0.1) / 0.1 : (abs($orig+0.1)) / 0.1;
		$BED->{$key}{odds} = myformat($BED->{$key}{odds});
		my $odds = $BED->{$key}{odds};
	}
	return $BED;
}

sub parse_fasta {
	my ($bed, $fastaFile, $outLog) = @_;
	open (my $in, "<", $fastaFile) or LOG ($outLog, "footPeak_kmer.pl: Cannot read from $fastaFile: $!\n") and return;
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
	my $GGG = 0;
	my $seq2 = $seq;
	while ($seq2 =~ /CG/g) {
		$CG ++;
	}
	undef($seq2);
	my $seq3 = $seq;
	while ($seq3 =~ /GGG/g) {
		$GGG ++;
	}
	undef($seq3);
	my $len = length($seq) - length($N);
	my $gc = int(100 * ($G+$C) / $len+0.5)/100;
	my $skew  = $gc == 0 ? 0 : int(100 * (abs($G-$C))  / ($G+$C)  +0.5)/100;
		$skew *= -1 if $C > $G;
	my $skew2 = $gc == 0 ? 0 : int(100 * (abs($G-$C)) / ($G+$C) + 0.5)/100;
	my $mod = $len == 0 ? 0 : $C > $G ? -1 * $C/$len : $G / $len;
		$skew2 *= $mod;
	my $cpg = ($C == 0 or $G == 0) ? 1 : int(100 * ($CG) / (($C) * ($G)) * $len + 0.5)/100;
	my $ggg = $G == 0 ? 1 : myformat(($len**2 * $GGG) / $G**3);
	$bed->{$number . 'a_len'}{$type} = $len;
	$bed->{$number . 'b_cpg'}{$type} = $cpg;
	$bed->{$number . 'c_gc'}{$type} = $gc;
	$bed->{$number . 'd_skew'}{$type} = $skew;
	$bed->{$number . 'e_skew2'}{$type} = $skew2;
	$bed->{$number . 'f_ggg'}{$type} = $ggg;
	$bed->{$number . 'h_acgt'}{$type} = "CG=$CG, C=$C, G=$G, A=$A, T=$T, len=$len";
	return $bed;
}


sub calculate_kmer { #TBA
	my ($seq) = @_;
}
