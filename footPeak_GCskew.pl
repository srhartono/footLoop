#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_n $opt_i);
getopts("vn:i:");

die "\nusage: $YW$0$N -i $LGN<geneIndexFile>$N -n $CY<footPeak output folder>$N\n\n" unless defined $opt_n and defined $opt_i and -e $opt_i and -d $opt_n;
my ($indexFile, $footPeak_folder) = ($opt_i, $opt_n);

my %gene;
#my ($indexFile) = "/group/stella/Work/Data/Fastq/110101_pacbio/1_Bed/180407_task1_GCskew/0_all_indexes_Alpha301.bed";
my @line = `cat $indexFile`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand, $feature) = split("\t", $line);
	die "Undefine gene at line=$line\n" if not defined $gene;
	$gene{$gene}{feature} = $feature;
}

makedir("$footPeak_folder/footPeak_GCskew/") if not -d "$footPeak_folder/footPeak_GCskew/";
makedir("$footPeak_folder/footPeak_GCskew/.TEMP") if not -d "$footPeak_folder/footPeak_GCskew/.TEMP";

open (my $outLog, ">", "$footPeak_folder/footPeak_GCskew/footPeak_GCskew_logFile.txt") or die "Failed to write to $footPeak_folder/footPeak_GCskew/footPeak_GCskew_logFile.txt: $!\n";

my ($footPeak_logFile) = "$footPeak_folder/footPeak_logFile.txt";
my $footLoop_run_script = `grep -iP "footLoop Run script\\s*:.+-g .+.fa" $footPeak_logFile`;
DIELOG($outLog, "Cannot find footLoop_run_script from footPeak logfile $footPeak_logFile\n") if not defined $footLoop_run_script or (defined $footLoop_run_script and $footLoop_run_script !~ /\w+/);
my @footLoop_run_script = split(" ", $footLoop_run_script);
my $genomeFile;
for (my $i = 1; $i < @footLoop_run_script; $i++) {
	next unless $footLoop_run_script[$i-1] eq "-g";
	$genomeFile = $footLoop_run_script[$i];
	last;
}
DIELOG($outLog, "Cannot find genome file from footPeak logfile $footPeak_logFile\n") if not defined $genomeFile or (defined $genomeFile and $genomeFile !~ /\w+/);
print $outLog "\ngenomeFile = $LCY$genomeFile$N\n\n";


# Get peaks from PEAKS_GENOME
my $cluster;
my @bedFiles = <$footPeak_folder/PEAKS_GENOME/*.genome.bed>;
foreach my $bedFile (sort @bedFiles) {
#	next if $bedFile !~ /CALM3_Pos.+CH/;
	$cluster = get_cluster($bedFile, $cluster, $outLog);
	preprocess_bed($bedFile, $outLog);
}
system("run_script_in_paralel2.pl \"fastaFromBed -fi $genomeFile -bed FILENAME -fo FILENAME.fa -s -name\" $footPeak_folder/footPeak_GCskew/.TEMP/ \"_[ABCDEFW].temp\" 1") == 0 or DIELOG($outLog, "Failed to run fastaFromBed: $!\n");
system("run_script_in_paralel2.pl \"rename.pl FILENAME PCB .PCB\" $footPeak_folder/footPeak_GCskew/.TEMP/ temp 1") == 0 or DIELOG($outLog, "Failed to run rename.pl: $!\n");
system("run_script_in_paralel2.pl \"counter_cpg_indiv.pl -w 200 -s 1 -o $footPeak_folder/footPeak_GCskew/ -A FILENAME\" $footPeak_folder/footPeak_GCskew/.TEMP/ _100.+temp.fa 1") == 0 or DIELOG($outLog, "Failed to run counter_cpg_indiv.pl: $!\n");
####### PARAMETERS
sub get_cluster {
	my ($bedFile, $cluster, $outLog) = @_; 
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($clusterFile) = "$footPeak_folder/FOOTCLUST/.TEMP/$bedFilename";
	$clusterFile =~ s/.genome.bed/.local.bed.clust/;
	if (-e $clusterFile) {
		my $linecount = -1;
		open (my $clusterIn, "<", $clusterFile) or DIELOG($outLog, "Failed to read from clusterFile $clusterFile: $!\n");
		while (my $line = <$clusterIn>) {
			chomp($line);
			$linecount ++;
			next if $linecount == 0;
			my ($id, $x, $xmax, $y, $ymax, $clust) = split("\t", $line);
			my ($id2, $number) = $id =~ /^(\d+)\.(\d+)$/;
			($id2) = $id if not defined $id2;
			print "id=$id, id2=$id2, number=$number, $line\n" if $linecount < 10;
			$id = $id2;
			if (not defined $cluster->{$id}) {
				$cluster->{$id}{clust} = $clust;
				$cluster->{$id}{len} = ($xmax - $x);
			}
			elsif (defined $cluster->{$id} and $cluster->{$id}{len} < $xmax - $x) {
				print "id=$id num=$number len=$cluster->{$id}{len} < $xmax-$x\n";
				$cluster->{$id}{clust} = $clust;
				$cluster->{$id}{len} = ($xmax - $x);
			}
		}
	}
	return($cluster);
}
sub preprocess_bed {
	my ($bedFile, $outLog) = @_; 
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my $window = 100;
	my $window2 = 200;
	my $outputA = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_A.temp";
	my $outputB = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_B.temp";
	my $outputC = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_C.temp";
	my $outputD = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_D.temp";
	my $outputE = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_E.temp";
	my $outputF = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_F.temp";
	my $outputW = "$footPeak_folder/footPeak_GCskew/.TEMP/$bedFilename\_$window\_W.temp";

	system("bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFile -o $outputA > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x -$window -y $window -i $bedFile -o $outputB > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x 0 -y $window2 -i $bedFile -o $outputC > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -x 0 -y 0 -i $bedFile -o $outputW > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window2 -y 0 -i $bedFile -o $outputD > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window -y $window -i $bedFile -o $outputE > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x 0 -y $window2 -i $bedFile -o $outputF > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");

#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window2 -y 0 -i FILENAME -o $outputA\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window -y $window -i FILENAME -o $outputB\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x 0 -y $window2 -i FILENAME -o $outputC\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
##	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -x 0 -y 0 -i FILENAME -o $outputW\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window2 -y 0 -i FILENAME -o $outputD\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window -y $window -i FILENAME -o $outputE\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x 0 -y $window2 -i FILENAME -o $outputF\" . \"(CALM3_Pos.+CH|FUS_Pos.+CH|MRFAP1L1_Neg.+GH).PEAK.genome.bed\" 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window2 -y 0 -i FILENAME -o $outputA\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window -y $window -i FILENAME -o $outputB\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x 0 -y $window2 -i FILENAME -o $outputC\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -x 0 -y 0 -i FILENAME -o $outputW\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window2 -y 0 -i FILENAME -o $outputD\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window -y $window -i FILENAME -o $outputE\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	system("run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x 0 -y $window2 -i FILENAME -o $outputF\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");

#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window2 -y 0 -i FILENAME -o $outputA\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x -$window -y $window -i FILENAME -o $outputB\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -a -x 0 -y $window2 -i FILENAME -o $outputC\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -x 0 -y 0 -i FILENAME -o $outputW\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window2 -y 0 -i FILENAME -o $outputD\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x -$window -y $window -i FILENAME -o $outputE\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
#	print "run_script_in_paralel2.pl \"bedtools_bed_change.pl -b -x 0 -y $window2 -i FILENAME -o $outputF\" $footPeak_folder/PEAKS_GENOME/ PEAK.genome.bed 1\n";# or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");

}


my @files = <$footPeak_folder/footPeak_GCskew/*.tsv>;
my %data;
my @header = ("label", "gene", "strand", "window", "threshold", "convtype", "wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$footPeak_folder/footPeak_GCskew/$N!\n" and exit if (@files == 0);
foreach my $input1 (sort @files) {
	my ($WINDOW, $SAMPLE, $TYPE);
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");
	#my ($label, $barcode, $desc, $gene, $strand, $window, $threshold, $convtype, $wind2, $sample, $type) = $fileName1 =~ /^(PCB.+)_(BC\d+)?_?(\w+)?_?gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	my @arr = $fileName1 =~ /^(PCB.+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	$arr[0] =~ s/^(PCB\d+)_.+$/$1/;
	my $outName = join("_", @arr[0..6]) . "_" . $arr[8];
	for (my $i = 0; $i < @arr; $i++) {
		die "Undefined i=$i header=$header[$i] arr[i] undef\n" if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		$data{data}{$outName}{$header[$i]} = $arr[$i];
		$WINDOW = $arr[$i] if $header[$i] eq "wind2";
		$SAMPLE = $arr[$i] if $header[$i] eq "sample";
		$TYPE = $arr[$i] if $header[$i] eq "type";
	}
	print "$input1:$LCY$outName$N\n";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	my $linecount = 0;
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		$linecount ++;
		my ($read, $value) = split("\t", $line);
		my ($id1, $id2, $id3) = $read =~ /^.*m(\d+_\d+)_\d+_c\w+_\w+_\w+\/(\d+)\/(ccs|\d+_\d+)/;
		$id3 = 0 if not defined $id3;
		$id3 = 0 if $id3 eq "ccs";
		DIELOG($outLog, "Failed to parse id1/2/3 from read=$LCY$read$N, file=$LGN$input1$N\n") if not defined $id1 or not defined $id2 or not defined $id3;
		my $id = "$id1$id2$id3"; $id =~ s/_//g;
		my $cluster = $cluster->{$id}{clust}; $cluster = -1 if not defined $cluster;
		print "id=$id, cluster=$cluster\n" if $linecount == 1;
		$data{id}{$outName}{$read} = $cluster;
		$data{read}{$outName}{$read}{$SAMPLE} = $value;
		$data{input}{$outName}{$SAMPLE} = join("_", @arr[0..5]);
	}	
	close $in1;
}

open (my $out1, ">", "$footPeak_folder/footPeak_GCskew/RESULT.TSV") or die "Cannot write to $footPeak_folder/footPeak_GCskew/RESULT.TSV: $!\n";
foreach my $outName (sort keys %{$data{input}}) {
	#open (my $out1, ">", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	print $out1 "file\tread\twindow\ttype\tfeature";
	foreach my $sample (sort keys %{$data{input}{$outName}}) {
		print $out1 "\t$sample";
	}
	print $out1 "\tcluster\n";
	last;
	#close $out1;
}

foreach my $outName (sort keys %{$data{read}}) {
#open (my $out1, ">>", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	my $GENE = $data{data}{$outName}{gene};
	my $feature = $gene{$GENE}{feature}; 
	$feature = "FEATURE_UNKNOWN" if not defined $feature;#print "Undef gene=$GENE feature\n" and next if not defined $feature;
	foreach my $read (sort keys %{$data{read}{$outName}}) {
		print $out1 "$data{input}{$outName}{$SAMPLE}\t$read\t$WINDOW\t$TYPE\t$feature";
		foreach my $sample (sort keys %{$data{read}{$outName}{$read}}) {
			print $out1 "\t$data{read}{$outName}{$read}{$sample}";
		}
		print $out1 "\t$data{id}{$outName}{$read}\n";
	}
	#close $out1;
}
__END__



close $out1;

__END__
PCB1_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed_100_E.temp.fa.dens.tsv
