#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_n $opt_i $opt_S $opt_G);
getopts("vn:i:SG:");

die "\nusage: $YW$0$N -i $LGN<geneIndexFile>$N -n $CY<footPeak output folder>$N\n\n" unless defined $opt_n and defined $opt_i and -e $opt_i and -d $opt_n;
my ($indexFile, $footPeak_folder) = ($opt_i, $opt_n);

my %gene;
my @line = `cat $indexFile`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand, $feature) = split("\t", $line);
	$feature = "FEATURE_UNKNOWN" if not defined $feature;
	die "Undefine gene at line=$line\n" if not defined $gene;
	$gene{$gene}{feature} = $feature;
}

my $outDir = "$footPeak_folder/footPeak_GCskew/";
makedir($outDir) if not -d $outDir;
makedir("$outDir/.TEMP") if not -d "$outDir.TEMP";

open (my $outLog, ">", "$outDir/footPeak_GCskew_logFile.txt") or die "Failed to write to $outDir/footPeak_GCskew_logFile.txt: $!\n";

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
my @files;
LOG($outLog, date() . "1. Getting cluster and preprocessing bed files!\n");
foreach my $bedFile (sort @bedFiles) {
	if (defined $opt_G) {next if $bedFile !~ /$opt_G/};
	my ($bedFileName) = $bedFile =~ /^$footPeak_folder\/PEAKS_GENOME\/(.+.PEAK.genome.bed)$/;
	my $tempFile;
	my @alphabets = qw(A B C W D E F);
	my @treats = qw(dens cont skew);
	for (my $i = 0; $i < @alphabets; $i++) {
		for (my $j = 0; $j < @treats; $j++) {
			$tempFile = "$footPeak_folder/footPeak_GCskew/$bedFileName";
			my $tsvFile = "$footPeak_folder/footPeak_GCskew/$bedFileName\_100\_$alphabets[$i].temp.fa.$treats[$j].tsv";
			#print "\t- $LCY$tsvFile$N\n";
			@files = (@files, $tsvFile);
		}
	}
	$tempFile =~ s/\/+/\//g;
	$cluster = get_cluster($bedFile, $tempFile, $cluster, $outLog);
	#print "$LCY$bedFile$N\n";
	preprocess_bed($bedFile, $outLog) if not defined $opt_S;
}
LOG($outLog, date() . "2. Calculating GC skew (might take a couple minutes)\n");
if (not defined ($opt_S)) {
	#print "run_script_in_paralel2.pl \"fastaFromBed -fi $genomeFile -bed FILENAME -fo FILENAME.fa -s -name\" $outDir.TEMP/ \"_[ABCDEFW].temp\" 1 > $outDir/.TEMP/fastaFromBed.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run fastaFromBed: $!\n");
	#print "run_script_in_paralel2.pl \"rename.pl FILENAME PCB .PCB\" $outDir.TEMP/ temp 1  > $outDir/.TEMP/rename.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run rename.pl: $!\n");
	#print "run_script_in_paralel2.pl \"counter_cpg_indiv.pl -w 200 -s 1 -o $outDir -A FILENAME\" $outDir.TEMP/ _100.+temp.fa 1  > $outDir/.TEMP/counter_cpg_indiv.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run counter_cpg_indiv.pl: $!\n");
	system("run_script_in_paralel2.pl \"fastaFromBed -fi $genomeFile -bed FILENAME -fo FILENAME.fa -s -name\" $outDir.TEMP/ \"_[ABCDEFW].temp\" 1 > $outDir/.TEMP/fastaFromBed.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run fastaFromBed: $!\n");
	system("run_script_in_paralel2.pl \"rename.pl FILENAME PCB .PCB\" $outDir.TEMP/ temp 1  > $outDir/.TEMP/rename.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run rename.pl: $!\n");
	system("run_script_in_paralel2.pl \"counter_cpg_indiv.pl -w 200 -s 1 -o $outDir -A FILENAME\" $outDir.TEMP/ _100.+temp.fa 1  > $outDir/.TEMP/counter_cpg_indiv.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run counter_cpg_indiv.pl: $!\n");
}
####### PARAMETERS
sub get_cluster {
	my ($bedFile, $tempFile, $cluster, $outLog) = @_; 
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($clusterFile) = "$footPeak_folder/FOOTCLUST/.TEMP/$bedFilename";
	$clusterFile =~ s/.genome.bed/.local.bed.clust/;
	my $maxClust = -1;
	if (-e $clusterFile) {
		my $linecount = -1;
		LOG($outLog, date() . "   bedFile=$LCY$bedFile$N, clusterFile=$LGN$clusterFile$N\n","NA");
		open (my $clusterIn, "<", $clusterFile) or DIELOG($outLog, "Failed to read from clusterFile $clusterFile: $!\n");
		while (my $line = <$clusterIn>) {
			chomp($line);
			$linecount ++;
			next if $linecount == 0;
			my ($id, $x, $xmax, $y, $ymax, $clust) = split("\t", $line);
			$maxClust = $clust if $maxClust < $clust;
			my ($id2, $number) = $id =~ /^(\d+)\.(\d+)$/;
			($id2) = $id if not defined $id2;
#			print "id=$id, id2=$id2, number=$number, $line\n" if $linecount < 10;
			$id = $id2;
			if (not defined $cluster->{$tempFile}{$id}) {
				my $currlen = $xmax - $x;
				#print "$clusterFile id=$id clust=$clust num=$number len=$currlen\n$tempFile,$id,$clust\n" if $id eq "1704132300151533640";
				$cluster->{$tempFile}{$id}{clust} = $clust;
				$cluster->{$tempFile}{$id}{len} = ($xmax - $x);
			}
			elsif (defined $cluster->{$tempFile}{$id} and $cluster->{$tempFile}{$id}{len} < $xmax - $x) {
				my $currlen = $xmax - $x;
				LOG($outLog, date() . " $clusterFile id=$id num=$number len=$cluster->{$tempFile}{$id}{len} < $currlen\n","NA");# if $id eq "1704132300151533640";
				$cluster->{$tempFile}{$id}{clust} = $clust;
				$cluster->{$tempFile}{$id}{len} = ($xmax - $x);
			}
		}
	}
	LOG($outLog, date() . "$LCY$bedFile$N cluster=$LGN$maxClust$N\n","NA");
	#print "temp=$LPR$tempFile$N\n";
	return($cluster);
}

sub preprocess_bed {
	my ($bedFile, $outLog) = @_; 
	LOG($outLog, date() . " Getting fasta and calculating GC skew from: $LCY$bedFile$N\n", "NA");
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my $window = 100;
	my $window2 = 200;
	my $outputA = "$outDir/.TEMP/$bedFilename\_$window\_A.temp";
	my $outputB = "$outDir/.TEMP/$bedFilename\_$window\_B.temp";
	my $outputC = "$outDir/.TEMP/$bedFilename\_$window\_C.temp";
	my $outputD = "$outDir/.TEMP/$bedFilename\_$window\_D.temp";
	my $outputE = "$outDir/.TEMP/$bedFilename\_$window\_E.temp";
	my $outputF = "$outDir/.TEMP/$bedFilename\_$window\_F.temp";
	my $outputW = "$outDir/.TEMP/$bedFilename\_$window\_W.temp";

=comment
	print "bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFile -o $outputA > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -a -x -$window -y $window -i $bedFile -o $outputB > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -a -x 0 -y $window2 -i $bedFile -o $outputC > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -x 0 -y 0 -i $bedFile -o $outputW > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -b -x -$window2 -y 0 -i $bedFile -o $outputD > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -b -x -$window -y $window -i $bedFile -o $outputE > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n\n";#;
	print "bedtools_bed_change.pl -b -x 0 -y $window2 -i $bedFile -o $outputF > $outputA.LOG 2>&1\n";# == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n";#;
=cut
	system("bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFile -o $outputA > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x -$window -y $window -i $bedFile -o $outputB > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x 0 -y $window2 -i $bedFile -o $outputC > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -x 0 -y 0 -i $bedFile -o $outputW > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window2 -y 0 -i $bedFile -o $outputD > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window -y $window -i $bedFile -o $outputE > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x 0 -y $window2 -i $bedFile -o $outputF > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
}

#my @files;# = <$outDir/PCB*.tsv>;
#if (defined $opt_G) {
#	@files = <$outDir/PCB*$opt_G*.tsv>;
#}
#else {
#	@files = <$outDir/PCB*.tsv>;
#}
LOG($outLog, date() . "Processing " . scalar(@files) . " files in $LCY$outDir$N\n");
my %data;
my @header = ("label", "gene", "strand", "window", "threshold", "convtype", "wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$outDir/$N!\n" and exit if (@files == 0);
foreach my $input1 (sort @files) {
	my ($tempFile) = $input1 =~ /^(.+)_100_.\.temp.fa.\w+.tsv$/;
	$tempFile =~ s/\/+/\//g;
	#print "input1=$input1, tempFile=$LCY$tempFile$N\n";
	$input1 =~ s/\/+/\//g;
	my ($WINDOW, $SAMPLE, $TYPE);
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");
	if (defined $opt_G) {
		next if $fileName1 !~ /$opt_G/;
	}
#	print "$input1\n";
	#my ($label, $barcode, $desc, $gene, $strand, $window, $threshold, $convtype, $wind2, $sample, $type) = $fileName1 =~ /^(PCB.+)_(BC\d+)?_?(\w+)?_?gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	next if $fileName1 !~ /^PCB/;
	my @arr = $fileName1 =~ /^(PCB.+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	$arr[0] =~ s/^(PCB\d+)_.+$/$1/;
	if (not defined $arr[0]) {
		for (my $i = 0; $i < @arr; $i++) {
			DIELOG($outLog, date() . "fileName=$fileName1 Undefined i=$i header=$header[$i] arr[i] undef\n") if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		}
	}
	my $outName = join("_", @arr[0..6]) . "_" . $arr[8];
	for (my $i = 0; $i < @arr; $i++) {
		DIELOG($outLog, date() . "Undefined i=$i header=$header[$i] arr[i] undef\n") if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		$data{data}{$outName}{$header[$i]} = $arr[$i];
		$WINDOW = $arr[$i] if $header[$i] eq "wind2";
		$SAMPLE = $arr[$i] if $header[$i] eq "sample";
		$TYPE = $arr[$i] if $header[$i] eq "type";
	}
	LOG($outLog, date() . "  - gene=$LGN$arr[1]$N pos=$LGN$arr[7]$N type=$LPR$arr[8]$N input=$outName$N\n","NA") if $arr[8] eq "skew";
	open (my $in1, "<", $input1) or DIELOG($outLog, date() . " Cannot read from $input1: $!\n");
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
		my $cluster = $cluster->{$tempFile}{$id}{clust}; $cluster = -1 if not defined $cluster;
		#print "id=$id, cluster=$cluster\n$input1,$id,$cluster\n" if $id eq "1704132300151533640";
		$data{cluster}{$outName}{$read} = $cluster;
		$data{id}{$outName}{$read} = $id;
		$data{read}{$outName}{$read}{$SAMPLE} = $value;
		$data{input}{$outName}{$SAMPLE} = join("_", @arr[0..5]);
	}	
	close $in1;
}
open (my $out1, ">", "$outDir/RESULT.TSV") or DIELOG($outLog, date() . "Cannot write to $outDir/RESULT.TSV: $!\n");
foreach my $outName (sort keys %{$data{input}}) {
	#open (my $out1, ">", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	print $out1 "file\tid\tcluster\tread\twindow\ttype\tfeature";
	foreach my $sample (sort keys %{$data{input}{$outName}}) {
		print $out1 "\t$sample";
	}
	print $out1 "\n";
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
		print $out1 "$data{input}{$outName}{$SAMPLE}\t$data{id}{$outName}{$read}\t$data{cluster}{$outName}{$read}\t$read\t$WINDOW\t$TYPE\t$feature";
		foreach my $sample (sort keys %{$data{read}{$outName}{$read}}) {
			print $out1 "\t$data{read}{$outName}{$read}{$sample}";
		}
		print $out1 "\n";
	}
	#close $out1;
}

GCskew_Rscript($outDir, $outLog);

sub GCskew_Rscript {
	my ($outDir, $outLog) = @_;
my $outDirFullpath = getFullpath($outDir);

my $RESULT = $outDir . "/RESULT.TSV";
my $LABEL = `cat $outDir/../.LABEL`; DIELOG($outLog, "Cannot find $outDir/../.LABEL!\n") if not defined $LABEL; chomp($LABEL);
$LABEL = $outDirFullpath . "/$LABEL";
my $Rscript = "
.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\",
\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\",
.libPaths()) )
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)    
library(ggplot2);
library(reshape2);

RESULT=\"$RESULT\";
PDFCLUST = paste(\"$LABEL\_BYCLUST.pdf\",sep=\"\")
PDFGENES = paste(\"$LABEL\_BYGENES.pdf\",sep=\"\")

df = read.table(RESULT,header=T,sep=\"\t\")
dm = melt(df,id.vars=c(\"file\",\"id\",\"cluster\",\"read\",\"window\",\"type\",\"feature\"));
dm\$feature = factor(dm\$feature,levels=c(\"PROMOTER\",\"GENEBODY\",\"TERMINAL\",\"FEATURE_UNKNOWN\"));
dm\$variable = factor(dm\$variable,levels=c(\"A\",\"B\",\"C\",\"W\",\"D\",\"E\",\"F\"))
dm\$gene = paste(dm\$file,dm\$feature)
genes = unique(dm\$gene)
dm\$group = paste(dm\$cluster)
pdf(PDFCLUST,width=7,height=7)
for (i in 1:length(genes)) {
   temp = dm[dm\$gene == genes[i] & dm\$type == \"skew\",]
   clusterz = unique(temp\$cluster)
   clusterz=clusterz[order(clusterz)]
	print(clusterz);
   for (j in 1:length(clusterz)) {
      temp2 = temp[temp\$cluster == clusterz[j],]
		counts = dim(temp2)[1] / length(unique(temp2\$variable))
      p = ggplot(temp2,aes(variable,value)) + geom_boxplot(aes(fill=variable),outlier.shape=NA) +
      theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(-1,1)) +
      ylab(\"GC Skew\") + xlab(\"Sample\") + ggtitle(paste(genes[i], \"cluster\",clusterz[j],\"\\ntotal read:\",counts))
      print(p)
   }
}
dev.off()

pdf(PDFGENES,width=7,height=7)
for (i in 1:length(genes)) {
   temp = dm[dm\$gene == genes[i] & dm\$type == \"skew\",]
	counts = dim(temp)[1] / length(unique(temp\$variable))
   p = ggplot(temp,aes(variable,value)) + geom_boxplot(aes(fill=variable),outlier.shape=NA) +
   theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(-1,1)) +
   ylab(\"GC Skew\") + xlab(\"Sample\") + ggtitle(paste(genes[i],\"\\ntotal read:\",counts))
   print(p)
}
dev.off()
";

LOG($outLog, date() . " $YW Running R script $LCY$outDir/RESULT.R$N\n");
open (my $outR, ">", "$outDir/RESULT.R") or DIELOG($outLog, date() . " Failed to write to $LCY$outDir/RESULT.R$N: $!\n");
print $outR $Rscript;
close $outR;
system("run_Rscript.pl $outDir/RESULT.R > $outDir/RESULT.R.LOG 2>&1") == 0 or DIELOG($outLog, date() . " Failed to run $LCY$outDir/RESULT.R$N: $!\n");
my $tailR = `tail $outDir/RESULT.R.LOG`;

LOG($outLog, "\n\n" . date() . " ${LGN}SUCCESS on running $LCY$outDir/RESULT.R$YW.\nLast 5 rows of log message:$N\n$N$tailR


${YW}To Run R script manually, do:$N
run_Rscript.pl $outDir/RESULT.R


${YW}Outputs:$N

- $LGN#grouped by each gene:$N
$LCY$LABEL\_BYGENES.pdf$N

- $LGN#grouped by each gene and each cluster:$N
$LCY$LABEL\_BYCLUST.pdf$N


");
}

__END__



close $out1;

__END__
PCB1_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed_100_E.temp.fa.dens.tsv
