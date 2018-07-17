#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_n $opt_i $opt_S $opt_G);
getopts("vn:i:SG:");

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

#my $OPTS = "vp:"; getopts($OPTS);
#use vars   qw($opt_v $opt_p);
#my @VALS =   ($opt_v,$opt_p);
#my $MAINLOG = myFootLog::MAINLOG($0, \@VALS, $OPTS);

my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
my @version = `cd $footLoopScriptsFolder && git log | head `;
my $version = "UNKNOWN";
foreach my $line (@version[0..@version-1]) {
   if ($line =~ /^\s+V\d+\.?\d*\w*\s*/) {
      ($version) = $line =~ /^\s+(V\d+\.?\d*\w*)\s*/;
   }
}
if (not defined $version or (defined $version and $version eq "UNKNOWN")) {
   ($version) = `cd $footLoopScriptsFolder && git log | head -n 1`;
}
if (defined $opt_v) {
   print "\n\n$YW$0 $LGN$version$N\n\n";
   exit;
}

################
# ARGV Parsing #
###############

die "\nusage: $YW$0$N -i $LGN<geneIndexFile.feature>$N -n $CY<footPeak output folder>$N

${LGN}Options:$N
-S: Skip the first 2 steps (e.g. if you've ran it before and would like to reuse the same files, as these usually don't change)
-G <gene>: Only run files with this gene in the name

" unless defined $opt_n and defined $opt_i and -e $opt_i and -d $opt_n;
my ($indexFile, $footPeakFolder) = ($opt_i, $opt_n);
my ($genewant) = $opt_G if defined $opt_G;

# sanity check -n footPeakFolder

($footPeakFolder) = getFullpath($footPeakFolder);
my $footClustFolder = "$footPeakFolder/FOOTCLUST/";
my $footKmerFolder  = "$footPeakFolder/KMER/";
my $uuid = getuuid();
my ($user) = $homedir =~ /home\/(\w+)/;
my $date = date();


############
# LOG FILE #
############
open (my $outLog, ">", "$footPeakFolder/footPeak_GCprofile_logFile.txt") or die date() . ": Failed to create outLog file $footPeakFolder/footPeak_GCprofile_logFile.txt: $!\n";
#open (my $outLog, ">", "$outDir/footPeak_GCprofile_logFile.txt") or die "Failed to write to $outDir/footPeak_GCprofile_logFile.txt: $!\n";
LOG($outLog, ">$0 version $version\n");
LOG($outLog, ">UUID: $uuid\n", "NA");
LOG($outLog, ">Date: $date\n", "NA");
LOG($outLog, ">Run script: $0 -i $opt_i -n $opt_n\n", "NA");


##########
# OUTDIR #
##########
my $outDir = "$footPeakFolder/GCPROFILE/";
if (-e "$footPeakFolder/footPeak_GCskew/" and not -d $outDir) {
	LOG($outLog, "Renaming $footPeakFolder/footPeak_GCskew into $outDir\n");
	system("/bin/mv $footPeakFolder/footPeak_GCskew $outDir") == 0 or DIELOG($outLog, date() . " Failed to rename $footPeakFolder/footPeak_GCskew into $outDir: $!\n");
}
if (not -d $outDir) {
	makedir($outDir);
	makedir("$outDir/.TEMP") if not -d "$outDir.TEMP";
}

##############
# INDEX FILE #
##############

my %gene;
my @line = `cat $indexFile`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand, $feature) = split("\t", $line);
	$feature = "FEATURE_UNKNOWN" if not defined $feature;
	die "Undefine gene at line=$line\n" if not defined $gene;
	$gene{$gene}{feature} = $feature;
}

##########################
# PARSE FOOTPEAK LOGFILE #
##########################

my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
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


###############################
# Get peaks from PEAKS_GENOME #
###############################
my $cluster;
my @bedFiles = <$footPeakFolder/PEAKS_GENOME/*.genome.bed>;
my @files;
LOG($outLog, date() . "1. Getting cluster and preprocessing bed files!\n");
foreach my $bedFile (sort @bedFiles) {
	if (defined $opt_G) {next if $bedFile !~ /$opt_G/};
	my ($bedFileName) = $bedFile =~ /^$footPeakFolder\/PEAKS_GENOME\/(.+.PEAK.genome.bed)$/;
	my $tempFile;
	my @alphabets = qw(A B C W D E F);
	my @treats = qw(dens cont skew);
	for (my $i = 0; $i < @alphabets; $i++) {
		for (my $j = 0; $j < @treats; $j++) {
			$tempFile = "$footPeakFolder/GCPROFILE/$bedFileName";
			my $tsvFile = "$footPeakFolder/GCPROFILE/$bedFileName\_100\_$alphabets[$i].temp.fa.$treats[$j].tsv";
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
	system("run_script_in_paralel2.pl \"fastaFromBed -fi $genomeFile -bed FILENAME -fo FILENAME.fa -s -name\" $outDir.TEMP/ \"_[ABCDEFW].temp\" 1 > $outDir/.TEMP/fastaFromBed.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run fastaFromBed: $!\n");
	system("run_script_in_paralel2.pl \"rename.pl FILENAME PCB .PCB\" $outDir.TEMP/ temp 1  > $outDir/.TEMP/rename.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run rename.pl: $!\n");
	system("run_script_in_paralel2.pl \"counter_cpg_indiv.pl -w 200 -s 1 -o $outDir -A FILENAME\" $outDir.TEMP/ _100.+temp.fa 1  > $outDir/.TEMP/counter_cpg_indiv.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run counter_cpg_indiv.pl: $!\n");
}
####### PARAMETERS
sub get_cluster {
	my ($bedFile, $tempFile, $cluster, $outLog) = @_; 
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($clusterFile) = "$footPeakFolder/FOOTCLUST/.TEMP/$bedFilename";
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
LOG($outLog, date() . "3. Processing " . scalar(@files) . " files in $LCY$outDir$N\n");
my %data;
my @header = ("label", "gene", "strand", "window", "threshold", "convtype", "wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$outDir/$N!\n" and exit if (@files == 0);
foreach my $input1 (sort @files) {
	my ($tempFile) = $input1 =~ /^(.+)_100_.\.temp.fa.\w+.tsv$/;
	$tempFile =~ s/\/+/\//g;
	#print "input1=$input1, tempFile=$LCY$tempFile$N\n";
	$input1 =~ s/\/+/\//g;
	my ($WINDOW, $SAMPLE, $TYPE);
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
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

GCprofile_Rscript($outDir, $outLog);

sub GCprofile_Rscript {
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
PDFCLUSTCpGdens = \"$LABEL\_BYCLUST_CpGdens.pdf\"
PDFGENESCpGdens = \"$LABEL\_BYGENES_CpGdens.pdf\"
PDFCLUSTGCcont = \"$LABEL\_BYCLUST_GCcont.pdf\"
PDFGENESGCcont = \"$LABEL\_BYGENES_GCcont.pdf\"
PDFCLUSTGCskew = \"$LABEL\_BYCLUST_GCskew.pdf\"
PDFGENESGCskew = \"$LABEL\_BYGENES_GCskew.pdf\"
PDFCLUST = c(PDFCLUSTCpGdens,PDFCLUSTGCcont,PDFCLUSTGCskew)
PDFGENES = c(PDFGENESCpGdens,PDFGENESGCcont,PDFGENESGCskew)

df = read.table(RESULT,header=T,sep=\"\\t\")
dm = melt(df,id.vars=c(\"file\",\"id\",\"cluster\",\"read\",\"window\",\"type\",\"feature\"));
dm\$feature = factor(dm\$feature,levels=c(\"PROMOTER\",\"GENEBODY\",\"TERMINAL\",\"FEATURE_UNKNOWN\"));
dm\$variable = factor(dm\$variable,levels=c(\"A\",\"B\",\"C\",\"W\",\"D\",\"E\",\"F\"))
dm\$gene = paste(dm\$file,dm\$feature)
genes = unique(dm\$gene)
dm\$group = paste(dm\$cluster)
##
clustTot = dim(aggregate(dm\$window,by=list(dm\$gene,dm\$cluster),sum))[1]
genesTot = dim(aggregate(dm\$window,by=list(dm\$gene),sum))[1]
clustCount = as.data.frame(plyr::count(dm,c(\"cluster\",\"gene\")));
clustCount\$freq = clustCount\$freq / (length(unique(dm\$variable)) * length(unique(dm\$type))); colnames(clustCount)[3] = \"clustGroup\"
genesCount = as.data.frame(aggregate(clustCount\$clustGroup,by=list(clustCount\$gene),sum));colnames(genesCount) = c(\"gene\",\"genesGroup\");
dm = merge(dm,clustCount,by=c(\"cluster\",\"gene\"),all=T)
dm = merge(dm,genesCount,by=c(\"gene\"),all=T)
dm\$clustGroup = paste(dm\$file,\" (\",dm\$feature,\") cluster \",dm\$cluster,\" (\",dm\$clustGroup,\" reads)\",sep=\"\")
dm\$genesGroup = paste(dm\$file,\" (\",dm\$feature,\") (\",dm\$genesGroup,\" reads)\",sep=\"\")
##
types = c(\"dens\",\"cont\",\"skew\")
ylimsMin=c(0,0,-1)
ylimsMax=c(1.2,1,1)
ylines=c(0.6,0.5,0)
ylabs = c(\"CpG Density\",\"GC Content\",\"GC Skew\")
for (i in 1:length(PDFCLUST)) {
   pdf(PDFCLUST[i],width=7,height=7*clustTot)
   temp = dm[dm\$type == types[i],]
   p = ggplot(temp,aes(variable,value)) +
      geom_boxplot(aes(fill=variable),outlier.shape=NA) +
      theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin[i],ylimsMax[i])) +
      annotate(geom=\"segment\",x=0,xend=8,y=ylines[i],yend=ylines[i],lty=2) +
      facet_grid(clustGroup~.) +
      ylab(ylabs[i]) + xlab(\"Samples\")
   print(p)
   dev.off()
}

for (i in 1:length(PDFGENES)) {
   pdf(PDFGENES[i],width=7,height=7*genesTot)
   temp = dm[dm\$type == types[i],]
   p = ggplot(temp,aes(variable,value)) +
      geom_boxplot(aes(fill=variable),outlier.shape=NA) +
      theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin[i],ylimsMax[i])) +
      annotate(geom=\"segment\",x=0,xend=8,y=ylines[i],yend=ylines[i],lty=2) +
      facet_grid(genesGroup~.) +
      ylab(ylabs[i]) + xlab(\"Samples\")
   print(p)
   dev.off()
}
";

LOG($outLog, date() . "4. $YW Running R script $LCY$outDir/RESULT.R$N\n");
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
$LCY$LABEL\_BYGENES_CpGdens.pdf$N
$LCY$LABEL\_BYGENES_GCcont.pdf$N
$LCY$LABEL\_BYGENES_GCskew.pdf$N

- $LGN#grouped by each gene and each cluster:$N
$LCY$LABEL\_BYCLUST_CpGdens.pdf$N
$LCY$LABEL\_BYCLUST_GCcont.pdf$N
$LCY$LABEL\_BYCLUST_GCskew.pdf$N


");
}

__END__



close $out1;

__END__
PCB1_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed_100_E.temp.fa.dens.tsv
