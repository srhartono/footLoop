#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_t $opt_R $opt_D);
getopts("vd:n:G:t:RD:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
}

use myFootLib;
use FAlite;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

my $date = getDate();
my $uuid = getuuid();

my ($dist, $footPeakFolder) = ($opt_d, $opt_n);
my $toggleRstrand = defined $opt_R ? "Yes" : "No";

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N $LGN-g gene$N $CY-n <footPeak's output folder (footPeak's -o)>$N

-D: tightness of read placement in each cluster (default: 200)
-R: toggle to order reads reversely (good foor in vitro REVERSE_ genes)

";

die $usage  unless defined $opt_n and -d $opt_n;

my $outDir = "$footPeakFolder/FOOTCLUST/";
$outDir= getFullpath($outDir);
my $windowDiv = defined $opt_D ? $opt_D : 200;
die "\nWindowDiv must be integer!\n" if $windowDiv !~ /^\d+$/ or $windowDiv <= 0;

my $clustThreshold = defined $opt_t ? $opt_t : 3;
die "\n\nERROR! -t must be positive number (pls don't make it scientific format, use a decimal)!\n" if ($clustThreshold !~ /^\d+\.?\d*$/);
makedir($outDir);
die "Failed to create output directory $LCY$outDir$N!\n" unless -d $outDir;
makedir("$outDir/.TEMP");
die "Failed to create output directory $LCY$outDir/.TEMP/$N!\n" unless -d "$outDir/.TEMP";
makedir("$outDir/PNG/");
die "Failed to create output directory $LCY$outDir/PNG$N!\n" unless -d "$outDir/PNG/";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footClust.txt") or die "Failed to create outLog file $outDir/logFile_footClust.txt: $!\n";


LOG($outLog, "

If R dies for any reason, make sure you have these required R libraries:
- RColorBrewer v1.1-2
- gridExtra v2.3
- labeling v0.3
- reshape2 v1.4.3
- ggplot2 v3.1.0
- GMD v0.3.3

");


# parse footLoop log file
($footPeakFolder) = getFullpath($footPeakFolder);
my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
my ($footLoopFolder);
my ($geneIndexFile);
open (my $infootPeak_logFile, "<", $footPeak_logFile) or DIELOG($outLog, "Failed to read from $footPeak_logFile: $!\n");
while (my $line = <$infootPeak_logFile>) {
	chomp($line);
		#>Options from footLoop.pl logfile=PCB9/0_Fastq/170804_pcb09_BSCset2_ccs_3minFP.fastq.gz.rmdup.fq.gz_MAP85pBUF100/logFile.txt	
	if ($line =~ /^Run script full/) {
		my @runscriptfull = split(" ", $line);
		for (my $i = 0; $i < @runscriptfull; $i++) {
			if ($runscriptfull[$i] eq "-n") {
				($footLoopFolder) = $runscriptfull[$i+1];
				last;

			}
		}
		last;
	}
}
my $footLoopFolderPrint = defined $footLoopFolder ? $footLoopFolder : "UNDEFINED";
die "Cannot find footLoop Folder $LCY$footLoopFolderPrint$N from logfile $footPeak_logFile\n" unless defined $footLoopFolder and -d $footLoopFolder;
close $infootPeak_logFile;
my $footLoop_logFile = $footLoopFolder . "/logFile.txt";
($geneIndexFile) = parse_footLoop_logFile($footLoop_logFile, $date, $uuid, $footLoopFolder);
my %coor;
open (my $inGeneIndexFile, "<", $geneIndexFile) or die;
while (my $line = <$inGeneIndexFile>) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	$coor{uc($gene)}{chr} = $chr;
	$coor{uc($gene)}{beg} = $beg;
	$coor{uc($gene)}{end} = $end;
	$coor{uc($gene)}{strand} = $strand;
}
close $inGeneIndexFile;
LOG($outLog, "cluster threshold = $clustThreshold\ngene Index file = $geneIndexFile\n");

# sanity check -d distance
$dist = 250 if not defined $opt_d;
LOG($outLog, date() . "\n$LRD Error$N: Distance must be digits!\n") and die unless $dist =~ /^\d+$/;

# get .fa file from footPeakFolder and copy
my ($faFile) = <$footPeakFolder/*.fa>;
LOG($outLog, date() . "\n$LRD Error$N: Cannot find any fasta file in $LCY$footPeakFolder$N\n\n") if not defined $faFile or not -e $faFile;
$faFile =~ s/\/+/\//g;
system("/bin/cp $faFile $outDir") == 0 or LOG($outLog, date() . "Failed to copy$LGN $faFile$N to $LCY$outDir$N: $!\n") and die;

# get .local.bed peak files used for clustering
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


# Parse fasta file
my %genes;
open (my $faIn2, "<", $faFile) or die;
open (my $faOut, ">", "$faFile.footClustfastafile") or die;
my $fasta2 = new FAlite($faIn2);
while (my $entry = $fasta2->nextEntry()) {
	my $def = uc($entry->def);
	my $seq = $entry->seq;
	print $faOut "$def\n$seq\n";
	my $length = length($seq);
	$def =~ s/>//;
	$genes{uc($def)} = $length;
}
close $faIn2;
close $faOut;
system("/bin/rm $faFile.footClustfastafile.fai") if -e "$faFile.footClustfastafile.fai";


####################
# Processing Input #
####################
my $label = "";
if (-e "$footPeakFolder/.LABEL") {
	($label) = `cat $footPeakFolder/.LABEL`;
	chomp($label);
}
else {
	DIELOG($outLog, "Failed to parse label from .LABEL in $footPeakFolder/.LABEL\n");
}
my $files_log = "#FOLDER=$footPeakFolder\n";
my $input1_count = -1;
my $curr_gene = -1;
my $type_count = -1;
foreach my $input1 (sort @local_peak_files) {
#DEBUG
	if (defined $opt_G and $input1 !~ /$opt_G/i) {
		LOG($outLog, date() . " Skipped $LCY$input1$N as it doesn't contain $LGN-G $opt_G$N\n");
		next;
	}

	# remove double // from folder
	$input1 =~ s/[\/]+/\//g;
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	my ($fullName1) = getFilename($input1, "full");

	# get gene and strand from file name
	my $parseName = parseName($fileName1);
   my ($label2, $gene, $strand, $window, $thres, $type) = @{$parseName->{array}};
   LOG($outLog, "\n------------- $YW WARNING!!! $N -------------\n\nUsing label=$label2. Inconsistent label in filename $LCY$fileName1$N\nLabel from $footPeakFolder/.LABEL: $label\nBut from fileName: $label2\n\n--------------- $YW WARNING!!! $N ---------------\n\n") if $label ne $label2;
   $label = $label2;
	my ($coorCHR, $coorBEG, $coorEND, $coorSTRAND) = ($coor{uc($gene)}{chr}, $coor{uc($gene)}{beg}, $coor{uc($gene)}{end}, $coor{uc($gene)}{strand});
	DIELOG($outLog, date() . "Fatal error at footClust.pl file=$input1, main gene index strand (\$coorSTRAND) is $coorSTRAND (has to be + or -\n") if $coorSTRAND !~ /^[\+\-]$/;
	my $Rstrand = $coorSTRAND eq "+" ? "Pos" : "Neg";#$coorSTRAND eq "-" ? "Neg" : "Pos";
#	my ($label2, $gene, $strand, $window, $thres, $type) = $fileName1 =~ /^(.+)\_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(GH|GC|CG|CH).PEAK/;
#	LOG($outLog, date() . "Cannot parse gene name from file=$LCY$input1$N\n") unless defined $gene and defined $strand and defined $window and defined $thres and defined $type;

	$thres *= 100 if $thres < 1;
	$gene   = uc($gene);

	if ($curr_gene ne $gene) {
		LOG($outLog, "\n$YW -------------- Processing $LGN$gene$N ------------- \n");
		$input1_count ++;
		$type_count = -1;
	}
	$type_count ++;
	$curr_gene = $gene;
	$strand = $strand eq "Neg" ? "-" : "+";
	# check peak file. Skip if there's less than 10 peaks
	my ($total_line) = `wc -l $input1` =~ /^(\d+)/;
	#if ($total_line <= 10) {
	#	LOG($outLog, "\n--------------\n\n$LGN$input1_count$N. ${LRD}NEXTED $N: $input1 ${LCY}because total peaks is less than 10 $N($LGN$total_line$N)\n\n");
	#	$files_log .= "$label\t$gene\t$strand\t$window\t$thres\t$type\t${LRD}Skip$N\ttotal_peak_all=$total_line\ttotal_read_unique=-1\ttotal_peak_used=-1\tPEAKS_LOCAL=PEAKS_LOCAL/$fullName1\tPEAK_FILE=NA\n";
	#	next;
	#}
	$files_log .= "$label\t$gene\t$strand\t$window\t$thres\t$type\t${LGN}Ran$N\ttotal_peak_all=$total_line";

	# Actual parse of peak file
	my ($linecount, $total_peak_used, $total_peak_all, %data) = (0,0,0);
	open (my $in1, "sort -k1,1 -k2,2n -k3,3n $input1|") or LOG($outLog, date() . "Cannot read from $input1: $!\n") and die;
	LOG($outLog, date() . "$LGN$input1_count.$type_count$N ${LCY}RUNNING$N: label=$LPR$label2$N gene=$LGN$gene$N strand=$LCY$strand$N window=$LGN$window$N thres=$LCY$thres$N type=$LGN$type$N input1=$input1\n");
	LOG($outLog, date() . "$LCY\tInfo$N: pcb=$label,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type\n");
	LOG($outLog, date() . "$LCY\tExample Lines$N:\n");
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my ($read, $beg, $end) = split("\t", $line);
		my ($id) = parse_readName($read, $outLog);#, $read =~ /^.*\.?m(\d+_\d+).+\/(\d+)\/(ccs|\d+_\d+)/;
		DIELOG($outLog, "\n\nERROR AT PARSING NUMBERS from read=$LGN$read$N\n\n") if not defined $id;

		LOG($outLog, date() . "$LRD\tERROR$N:$LCY Read must end in this format$N: <anything>/<hole number>/ccs\n\n$read\n\n") and die if not defined $id;
		my $check = 0;
		$total_peak_all ++;

		# reads with multiple peaks separated by less than distance (250) is merged together
		if (defined $data{$id}) {
			for (my $i = 0; $i < @{$data{$id}}; $i++) {
				my $beg2 = $data{$id}[$i][0];
				my $end2 = $data{$id}[$i][1];
				if ($beg < $end2 + $dist) {
					$data{$id}[$i][1] = $end;
					$check = 1;
					LOG($outLog, date() . "$LGN\tline=$linecount$N: readid=$id MERGED beg2=$beg2 end2=$end2 with beg=$beg end=$end into beg3=$beg2 end3=$end\n") if $linecount < 5;
					last;
				}
			}
		}

		# if 1 peak or peak is far then create new peak
		if ($check == 0) {
			push(@{$data{$id}}, [$beg,$end,$read,$gene]);
			$total_peak_used ++;
			LOG($outLog, date() . "$LGN\tline=$linecount$N: readid=$id beg=$beg end=$end\n") if $linecount < 5;
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
	my $tempFile1_Peak = "$outDir/.TEMP/.$fullName1.clust.temp1";
	my $tempFile2_Clust = "$outDir/.TEMP/.$fullName1.clust";
	my $tempFile2_PNG = "$outDir/PNG/$fullName1.clust.png";
	open (my $out1, ">", "$tempFile1_Peak") or DIELOG($outLog, date() . __LINE__ . "\tCannot write to $LCY$fileName1.temp$N: $!\n");
	print $out1 "id\txbeg\txend\ttotal_peak\n";
	foreach my $id (sort keys %data) {
		my $total_peak = @{$data{$id}}; #add
		for (my $i = 0; $i < @{$data{$id}}; $i++) {
			my $xbeg = $data{$id}[$i][0];
			my $xend = $data{$id}[$i][1];
			my $read = $data{$id}[$i][2];
			print $out1 "$id\t$xbeg\t$xend\t$total_peak\n";
		}
	}
	close $out1;

	# Write and run R script
	LOG($outLog, date() . "$LCY\tRunning R script$N $tempFile1_Peak.R\n");
	my $Rscript = "
print(.libPaths())
		require(stringi)
		library(labeling)
		library(ggplot2)
		library(reshape2)
		library(grid)
		library(gridExtra)
		library(RColorBrewer)
	library(Cairo)
	set.seed(420)
	df.main = read.table(\"$tempFile1_Peak\",header=T,sep=\"\\t\",colClasses=c(\"factor\",\"integer\",\"integer\",\"integer\")) # changed df to df.main
	df.main.total_peak = unique(df.main[,4])

	last_cluster = 0

	df = df.main
	df2 = df[,-1]
	lenz = length(unique(paste(df\$xbeg,df\$xend)))
	k.best = 0
	if (lenz >= 20) {
		set.seed(420)
		for (i in 2:20) {
			k = kmeans(df2,i)
			mylen = length(k[k\$withinss/1e6 > $clustThreshold])
			if (mylen == 0) {
				k.best = i
				break
			}
		}
		if (k.best == 0) {
			lenz = 5
		} else { lenz = k.best}
		df\$clust = k\$cluster
		print(paste(\"lenz=\",lenz,\"cluster center = \",k.best))
	} else if (lenz > 5) {
		k = kmeans(df2,2)
		df\$clust = k\$cluster
		print(paste(\"lenz=\",lenz,\"cluster center = 2\"))
	} else {
		df\$clust = 1
		print(paste(\"lenz=\",lenz,\"cluster center = 1\"))
	}


	Rstrand = \"\$Rstrand\"
	toggleRstrand = \"$toggleRstrand\"

	if (Rstrand == \"Neg\" & toggleRstrand == \"Yes\") {
		df = df[order(-df\$clust, -as.integer(df\$xend/$windowDiv), -as.integer(df\$xbeg/$windowDiv)),]
	} else {
		df = df[order(df\$clust, as.integer(df\$xbeg/$windowDiv), as.integer(df\$xend/$windowDiv)),]
	}

	# aggregate mean xbeg and xend position then sort cluster number
	df2 = as.data.frame(aggregate(subset(df,select=c(xbeg,xend)),by=list(df\$clust),function(x)mean(x,trim=0.05)))
	colnames(df2) = c(\"clust\",\"xbeg\",\"xend\")
	if (Rstrand == \"Neg\" & toggleRstrand == \"Yes\") {
		df2 = df2[order(-as.integer(df2\$xend/$windowDiv), -as.integer(df2\$xbeg/$windowDiv)),]
	} else {
		df2 = df2[order(as.integer(df2\$xbeg/$windowDiv), as.integer(df2\$xend/$windowDiv)),]
	}
	df2\$clust2 = seq(1,dim(df2)[1])

	# merge cluster number to df
	df = as.data.frame(merge(df,subset(df2,select=c(clust,clust2) ) ))
	df\$clust = df\$clust2
	df\$id2 = paste(df\$id,\".\",df\$clust,sep=\"\")
	df = subset(df,select=-clust2)

	mylen = $windowDiv;
	if (Rstrand == \"Neg\" & toggleRstrand == \"Yes\") {
		df = df[order(-df\$clust, -as.integer(df\$xend/mylen), -as.integer(df\$xbeg/mylen)),]
	} else {
		df = df[order(df\$clust, as.integer(df\$xbeg/mylen), as.integer(df\$xend/mylen)),]
	}
	df\$ybeg = seq(1,dim(df)[1])
	df\$yend = seq(2,dim(df)[1]+1)
	curr_last_cluster = max(df\$clust)
	df\$clust = df\$clust + last_cluster
	last_cluster = curr_last_cluster
	df = subset(df,select=c(id,xbeg,xend,ybeg,yend,clust,id2));# df\$h = h
	df.final = df

	df = df.final
	df\$ybeg = seq(1,dim(df)[1])
	df\$yend = seq(2,dim(df)[1]+1)
	png(type=\"cairo\",\"$tempFile2_PNG\",height=(dim(df)[1])*5,width=max(df\$xend)+10)
	ggplot(df, aes(xbeg,ybeg)) + geom_rect(aes(xmin=xbeg,ymin=ybeg,xmax=xend,ymax=yend,fill=as.factor(clust))) + theme_bw() + 
	theme(panel.grid=element_blank()) + coord_fixed(ratio=10)
	dev.off()
	write.table(df,\"$tempFile2_Clust\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
	";


	open (my $outR, ">", "$tempFile1_Peak.R") or DIELOG($outLog, date() . __LINE__ . "Failed to write to $tempFile1_Peak.R: $!\n");
	print $outR $Rscript;
	system("R --vanilla --no-save < $tempFile1_Peak.R > $tempFile1_Peak.R.log 2>&1");
	close $outR;
	my @Rlog = `tail -n 3 $tempFile1_Peak.R.log`;
	LOG($outLog, date() . "$LCY\tR logs$N:\n\t\t\t\t" . join("\t\t\t\t", @Rlog) . "\n");
	
	
	# print genome indiv

	my ($ampliconLength) = $genes{uc($gene)};
	DIELOG($outLog, date() . " Failed to get amplicon length frmo gene=" . uc($gene) . ": $!\n") if not defined $ampliconLength;

	my $outlocalindivBedFile = "$outDir/CLUST_LOCAL/$fullName1.indiv.clust.bed";
	my $outlocalindivClustFile = "$outDir/CLUST_LOCAL/$fullName1.indiv.clust";
	my $outgenomeindivBedFile = "$outDir/CLUST_GENOME/$fullName1.indiv.clust.bed";
	my $outgenomeindivClustFile = "$outDir/CLUST_GENOME/$fullName1.indiv.clust";
	$outgenomeindivBedFile =~ s/\.local/.genome/i;
	$outgenomeindivClustFile =~ s/\.local/.genome/i;
	my $outlocalBedFile  = "$outDir/CLUST_LOCAL/$fullName1.clust";
	my $outgenomeBedFile  = "$outDir/CLUST_GENOME/$fullName1.clust";
	$outgenomeBedFile =~ s/\.local/.genome/i;
	makedir("$outDir/CLUST_LOCAL") if not -d "$outDir/CLUST_LOCAL";
	makedir("$outDir/CLUST_GENOME") if not -d "$outDir/CLUST_GENOME";

	# process clust
	my %clust0;
	open (my $clustIn0, "<", $tempFile2_Clust) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $tempFile2_Clust: $!\n");
	while (my $line = <$clustIn0>) {
		chomp($line);
		if ($line =~ /xbeg\txend/) {next;}
		my ($id, $xbeg, $xend, $ybeg, $yend, $clust, $id2) = split("\t", $line);#, $total_peaks) = split("\t", $line);
		@{$clust0{$id}{clustarr}} = defined $clust0{$id} ? (@{$clust0{$id}{clustarr}},$clust) : ($clust);
		$clust0{$id}{total_peak} ++;
		push(@{$clust0{$id}{line}}, "$id\t$xbeg\t$xend");
		$clust0{$id}{beg} = (not defined $clust0{$id}{beg}) ? int($xbeg / $windowDiv) : $clust0{$id}{beg} < int($xbeg / $windowDiv) ? $clust0{$id}{beg} : int($xbeg / $windowDiv);
		$clust0{$id}{end} = (not defined $clust0{$id}{end}) ? int($xend / $windowDiv) : $clust0{$id}{end} > int($xend / $windowDiv) ? $clust0{$id}{end} : int($xend / $windowDiv);
		$clust0{$id}{id2} = $id2;
	}
	close $clustIn0;

	foreach my $id (sort {$clust0{$a}{total_peak} <=> $clust0{$b}{total_peak}} keys %clust0) {
		my $iter = 0;
		foreach my $h (sort {$a <=> $b} @{$clust0{$id}{clustarr}}) {
			$clust0{$id}{clust} += (100**$iter * $h);
			$iter ++;
		}
	}

#	if (Rstrand == \"Neg\" & toggleRstrand == \"Yes\") {
#		df = df[order(-df\$clust, -as.integer(df\$xend/mylen), -as.integer(df\$xbeg/mylen)),]
#	} else {
#		df = df[order(df\$clust, as.integer(df\$xbeg/mylen), as.integer(df\$xend/mylen)),]
#	}
	my ($ybeg, $yend) = (1,2);
	my $clust = 0; my $lastClustText = "NA";
	open (my $outClust0, ">", $tempFile2_Clust) or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $tempFile2_Clust: $!\n");

	if ($Rstrand eq "Neg" and $toggleRstrand eq "Yes") {
		foreach my $id (sort {$clust0{$b}{total_peak} <=> $clust0{$a}{total_peak} || $clust0{$b}{clust} <=> $clust0{$a}{clust} || $clust0{$b}{end} <=> $clust0{$a}{end} || $clust0{$b}{beg} <=> $clust0{$a}{beg} } keys %clust0) {
			my $clustText = join(",", @{$clust0{$id}{clustarr}});# . " (" . $clust0{$id}{total_peak} . ")";
			$clust ++ if $lastClustText ne $clustText;
			$lastClustText = $clustText;
		}
		$clust ++;
		$lastClustText = "NA";
		foreach my $id (sort {$clust0{$b}{total_peak} <=> $clust0{$a}{total_peak} || $clust0{$b}{clust} <=> $clust0{$a}{clust} || $clust0{$b}{end} <=> $clust0{$a}{end} || $clust0{$b}{beg} <=> $clust0{$a}{beg} } keys %clust0) {
			my $clustText = join(",", @{$clust0{$id}{clustarr}});# . " (" . $clust0{$id}{total_peak} . ")";
			$clust -- if $lastClustText ne $clustText;
			foreach my $line (@{$clust0{$id}{line}}[0..(@{$clust0{$id}{line}}-1)]) {
				print $outClust0 "$line\t$ybeg\t$yend\t$clust\t$clustText\n";
			}
			$ybeg ++;
			$yend ++;
			$lastClustText = $clustText;
		}
	}
	else {
		foreach my $id (sort {$clust0{$a}{total_peak} <=> $clust0{$b}{total_peak} || $clust0{$a}{clust} <=> $clust0{$b}{clust} || $clust0{$a}{end} <=> $clust0{$b}{end} || $clust0{$a}{beg} <=> $clust0{$b}{beg} } keys %clust0) {
			my $clustText = join(",", @{$clust0{$id}{clustarr}});# . " (" . $clust0{$id}{total_peak} . ")";
			$clust ++ if $lastClustText ne $clustText;
			foreach my $line (@{$clust0{$id}{line}}[0..(@{$clust0{$id}{line}}-1)]) {
				print $outClust0 "$line\t$ybeg\t$yend\t$clust\t$clustText\n";
			}
			$ybeg ++;
			$yend ++;
			$lastClustText = $clustText;
		}
	}
	close $outClust0;

	open (my $clustIn, "<", $tempFile2_Clust) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $tempFile2_Clust: $!\n");
	open (my $outlocalindivClust, ">", $outlocalindivClustFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outlocalindivClustFile: $!\n");
	open (my $outlocalindivBed, ">", $outlocalindivBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outlocalindivBedFile: $!\n");
	open (my $outgenomeindivClust, ">", $outgenomeindivClustFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outgenomeindivClustFile: $!\n");
	open (my $outgenomeindivBed, ">", $outgenomeindivBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outgenomeindivBedFile: $!\n");


	my %used; my %cl;
	LOG($outLog, date() . "$LCY\tCreating fasta file for each cluster$N $tempFile2_Clust.fa\n");
	while (my $line = <$clustIn>) {
		chomp($line);
		if ($line =~ /xbeg\txend/) {next;}
		my ($id, $xbeg, $xend, $ybeg, $yend, $clust, $id2) = split("\t", $line);
		die "Undefined clust at line:\n\n$line\n\n" if not defined $clust;
		my $xBEG = $coorBEG + $xbeg;
		my $xEND = $coorBEG + $xend;
		print $outlocalindivBed "$coorCHR\t$xbeg\t$xend\t$gene.$id\t$clust;$id2\t$coorSTRAND\n";
		print $outlocalindivClust "$id\t$xbeg\t$xend\t$ybeg\t$yend\t$clust\t$id2\n";
		print $outgenomeindivBed "$coorCHR\t$xBEG\t$xEND\t$gene.$id\t$clust;$id2\t$coorSTRAND\n";
		print $outgenomeindivClust "$id\t$xBEG\t$xEND\t$ybeg\t$yend\t$clust\t$id2\n";
		my ($xmid) = int(($xend + $xbeg)/2+0.5);
		push(@{$cl{$clust}{beg}}, $xbeg);
		push(@{$cl{$clust}{mid}}, $xmid);
		push(@{$cl{$clust}{end}}, $xend);
		$cl{$clust}{total_read} ++;
		$cl{$clust}{total_read_unique} ++ if not defined $used{$id};
		$used{$id} = 1;
	}
	#=commentcut2
	#open (my $outz2, ">", "$outDir/.TEMP/$fullName1.clust") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outDir/.TEMP/$fullName1.clust: $!\n");
	#print $outz2 "$clustheader\n";
	#foreach my $num (@clust[0..@clust-1]) {
	#	print $outz2 "$clust{$num}\n";
	#}
	#close $outz2;

	my $tempFile3_Bed = $tempFile2_Clust . ".bed";
	my $tempFile4_Fa  = $tempFile3_Bed . ".fa";
	open (my $out2, ">", "$tempFile3_Bed") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $tempFile3_Bed: $!\n");
	open (my $outlocalBed, ">", $outlocalBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outlocalBedFile: $!\n");
	open (my $outgenomeBed, ">", $outgenomeBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outgenomeBedFile: $!\n");
	print $out2 "#gene\tbeg\tend\tcluster\ttotal_peak\tstrand\n";
	#my ($coorCHR, $coorBEG, $coorEND, $coorSTRAND) = ($coor{uc($gene)}{chr}, $coor{uc($gene)}{beg}, $coor{uc($gene)}{end}, $coor{uc($gene)}{strand});
	DIELOG($outLog, date() . "\tFailed to get coordinate for gene=" . uc($gene) . ", chr=$coorCHR beg=$coorBEG end=$coorEND strand=$coorSTRAND\n") if not defined $coorCHR or not defined $coorBEG or not defined $coorEND or not defined $coorSTRAND;
	foreach my $clust (sort {$a <=> $b} keys %cl) {
		my $BEG   = int(tmm(@{$cl{$clust}{beg}})+0.5);
		my $MID   = int(tmm(@{$cl{$clust}{mid}})+0.5);
		my $END   = int(tmm(@{$cl{$clust}{end}})+0.5);
		my $BEGSD = int(tmmsd(@{$cl{$clust}{beg}})+0.5);
		my $MIDSD = int(tmmsd(@{$cl{$clust}{mid}})+0.5);
		my $ENDSD = int(tmmsd(@{$cl{$clust}{end}})+0.5);
		my ($begOfBEG, $endOfBEG) = ($BEG - $BEGSD, $BEG + $BEGSD);
		my ($begOfMID, $endOfMID) = ($MID - $MIDSD, $MID + $MIDSD);
		my ($begOfEND, $endOfEND) = ($END - $ENDSD, $END + $ENDSD);
		my ($begOfALL, $endOfALL) = ($begOfBEG, $endOfEND);
		if ($endOfEND > $ampliconLength) {
			LOG($outLog, date() . "    INFO: clust=$clust $gene BEG=$BEG END=$END begOfBEG=$begOfBEG endOfEND = $endOfEND changed to ampliconLength=$ampliconLength\n");
		}
		($begOfBEG) = 1 if $begOfBEG < 1;
		($endOfBEG) = $ampliconLength if $endOfBEG > $ampliconLength;
		($begOfMID) = 1 if $begOfMID < 1;
		($endOfMID) = $ampliconLength if $endOfMID > $ampliconLength;
		($begOfEND) = 1 if $begOfEND < 1;
		($endOfEND) = $ampliconLength if $endOfEND > $ampliconLength;
		($begOfALL) = 1 if $begOfALL < 1;
		($endOfALL) = $ampliconLength if $endOfALL > $ampliconLength;
#		my $total = @{$cl{$clust}{beg}};
		my $total_read = $cl{$clust}{total_read};
		my $total_read_unique = $cl{$clust}{total_read_unique};
		my $begGENOME = $coorBEG + $begOfALL;
		my $endGENOME = $coorBEG + $endOfALL;
		print $outgenomeBed "$coorCHR\t$begGENOME\t$endGENOME\t$gene.$clust\t$total_read_unique.$total_read\t$coorSTRAND\n";
		print $outlocalBed "$coorCHR\t$begOfALL\t$endOfALL\t$gene.$clust\t$total_read_unique.$total_read\t$coorSTRAND\n";
		print $out2 "$gene\t$begOfALL\t$endOfALL\t$clust.WHOLE\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfBEG\t$endOfBEG\t$clust.BEG\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfMID\t$endOfMID\t$clust.MID\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfEND\t$endOfEND\t$clust.END\t$total_read_unique.$total_read\t$strand\n";
	}
	close $out2;
	LOG($outLog, "$YW ::: fastaFromBed -fi $faFile.footClustfastafile -bed $tempFile3_Bed -fo $tempFile4_Fa -s ::: $N\n","NA");
	system("fastaFromBed -fi $faFile.footClustfastafile -bed $tempFile3_Bed -fo $tempFile4_Fa -s > $tempFile3_Bed.LOG 2>&1") == 0 or LOG($outLog, date() . " Failed to fastaFromBed -fi$LCY $faFile.footClustfastafile$N -bed$LGN $tempFile3_Bed$N -fo$LPR $tempFile4_Fa$N -s >$YW $tempFile3_Bed.LOG$N 2>\&1: $!\n");

	LOG($outLog, date() . "$LCY\tAnalyzing Sequence Content from$N $tempFile4_Fa\n");
	open (my $faIn, "$tempFile4_Fa") or DIELOG($outLog, "Failed to read from $tempFile4_Fa: $!\n");
	my $fasta = new FAlite($faIn); $linecount = 0;
	while (my $entry = $fasta->nextEntry()) {
		$linecount ++;
		my $def = uc($entry->def); $def =~ s/^>//;
		my $seq = uc($entry->seq());
		my ($seqshort) = $seq =~ /^(.{0,10})/; $seqshort = "" if not defined $seqshort;
		LOG($outLog, "def=$LCY>$def$N seq=$LGN$seqshort$N\n\n","NA") if $linecount == 1;
	}
}	
open (my $outz, ">", "$outDir/.0_LOG_FILESRAN") or DIELOG($outLog, "Failed tow rite to $outDir/.0_LOG_FILESRAN:$!\n");
print $outz $files_log;
close $outz;
LOG($outLog, "\n$YW ----------- FILES RAN ---------- $N\n\n" . date() . "$LCY\tSummary of run$N: $outDir/.0_LOG_FILESRAN\n\n$files_log\n\n");


sub parse_footLoop_logFile {
   my ($logFile, $date, $uuid, $footFolder, $version) = @_;
   #my @line = `cat $logFile`;
   my $paramsFile = "$footFolder/.PARAMS";
 #  my $inputFolder = $defOpts->{n};
	my $geneIndexFile;
   my @parline = `cat $paramsFile`;
   foreach my $parline (@parline) {
      if ($parline =~ /footLoop.pl,geneIndexFile,/) {
         ($geneIndexFile) = $parline =~ /geneIndexFile,(.+)$/;
      }
	}
	return($geneIndexFile);
}

__END__
=comment
      if ($parline =~ /footLoop.pl,geneIndexFile,/) {
			($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
      }
      if ($parline =~ /footLoop.pl,geneIndexFile,/) {
         ($defOpts->{geneIndexFile}) = $parline =~ /geneIndexFile,(.+)$/;
      }
      if ($parline =~ /footLoop.pl,seqFile,/) {
         ($defOpts->{seqFile}) = $parline =~ /seqFile,(.+)$/;
      }
      if ($parline =~ /footLoop.pl,samFile,/) {
         ($defOpts->{samFile}) = $parline =~ /samFile,(.+)$/;
      }
      if ($parline =~ /footLoop.pl,samFixedFile,/) {
         ($defOpts->{samFixed}) = $parline =~ /samFixedFile,(.+)$/;
      }
      if ($parline =~ /footLoop.pl,samFixedFileMD5,/) {
         ($defOpts->{samFixedMD5}) = $parline =~ /samFixedFileMD5,(.+)$/;
      }
   }
   # %others contain other parameters from footLoop that aren't options (e.g. uuid)
   my ($other, $outLog);


   foreach my $line (@line[0..@line-1]) {
      if ($line =~ /^\s*Date\s*:/) {
         ($other->{date}) = $line =~ /^Date\s+:\s*([a-zA-Z0-9\:\-].+)$/;
         print "Undefined date from $line\n" and die if not defined $other->{date};
      }
      if ($line =~ /^\s*Run ID\s*:/) {
         ($other->{uuid}) = $line =~ /^Run ID[ \t]+:[ \t]+([a-zA-Z0-9]+.+)$/;
         print "Undefined uuid from $line\n" and die if not defined $other->{uuid};
      }
      if ($line =~ /^\s*Run script\s*:/) {
         ($other->{footLoop_runscript}) = $line =~ /^Run script[ ]+:(.+)$/;
         print "Undefined runscript from $line\n" and die if not defined $other->{footLoop_runscript};
         ($defOpts, $other->{runscript}) = parse_runscript($defOpts, $usrOpts, $other->{footLoop_runscript});
      }
      if ($line =~ /^\s*Output\s*:/) {
         my ($value) = $line =~ /Output.+(\.0_orig_\w{32})/;
         die "Died at line=$line, value=?\n" if not defined $value;
         $value = $footFolder . "/$value";
         die "Died at line=$line, value=?\n" if not -d $value;
         ($other->{origDir}) = $value;
         print "Undefined origDir from input=$inputFolder, line=$line\n" and die if not defined $other->{origDir};
         print "origDir $other->{origDir} does not exist!\n" and die if not -d $other->{origDir};
         $defOpts->{origDir} = $other->{origDir};
         ($other->{md5}) = $other->{origDir} =~ /\.0_orig_(\w{32})/;
         print "Can't parse md5 from outdir (outdir=$defOpts->{origDir})\n" and die if not defined $other->{md5};
      }
#     if ($line =~ /geneIndexFile=/) {
#        ($defOpts->{geneIndexFile}) = $line =~ /geneIndexFile=(.+)$/ if $line !~ /,gene=.+,beg=\d+,end=\d+$/;
#        ($defOpts->{geneIndexFile}) = $line =~ /geneIndexFile=(.+),gene=.+,beg=\d+,end=\d+$/ if $line =~ /,gene=.+,beg=\d+,end=\d+$/;
#        $defOpts->{geneIndexFile} = $footFolder . "/" .  getFilename($defOpts->{geneIndexFile});
#     }

      if ($line =~ /^!\w+=/) {
         my ($param, $value) = $line =~ /^!(\w+)=(.+)$/;
         my $param2 = defined $param ? $param : "__UNDEF__";
         my $value2 = defined $value ? $value : "__UNDEF__";
         if ($value =~ /\//) {
            if ($value =~ /\/?\.0_orig\w{32}/) {
               ($value) = $value =~ /^.+(\/?\.0_orig_\w{32})/;
               $value = $footFolder . "/$value";
               die "Died at line=line, param=$param, value=?\n" if not defined $value;
            }
            if ($value =~ /\/\.geneIndex/) {
               ($value) = $value =~ /^.+(\/\.geneIndex.+)$/;
               $value = $footFolder . "/$value";
               die "Died at line=line, param=$param, value=?\n" if not defined $value;
            }
            else {
               ($value) = getFilename($value, 'full');
               $value = $footFolder . "/$value";
               die "Died at line=line, param=$param, value=?\n" if not defined $value;
            }
         }
         print "$param = $value\n";
         print "Cannot parse param=$param2 and value=$value2 from line=$line\n" and die if not defined $param or not defined $value;
         print "$param file $value does not exist!\n" and die if $value =~ /\/+/ and not -e $value;
         ($defOpts->{$param}) = $value if $param ne "n";
      }
   }
   $defOpts->{o} = $defOpts->{n} if not defined $opt_o;
   makedir($defOpts->{o}) if not -d $defOpts->{o};
   #die "opt = $opt_o = $defOpts->{o}\n";
   open ($outLog, ">", "$defOpts->{o}/footPeak_logFile.txt") or print "Failed to write to $defOpts->{o}/footPeak_logFile.txt: $!\n" and die;
   record_options($defOpts, $usrOpts, $usrOpts2, $other, $outLog, $logFile, $date, $uuid, $version);
#  print "Output = $defOpts->{o}\n";
   return($defOpts, $outLog);
}

=cut
#commentcut2
	my %cl;
	my %clust; my %clustlen; my @clust; my $clustheader = "";
	open (my $in2, "<", "$tempFile2_Clust") or DIELOG($outLog, date() . __LINE__ . "\tFailed to read from tempFile2_Clust: $!\n");
	my $linecount5 = 0;
	my %used;
	while (my $line = <$in2>) {
		$linecount5 ++;
		chomp($line);
		if ($line =~ /(clust|xmax|ymax)/) {$clustheader = $line;next;}
		my ($orignum, $beg, $end, $y, $y2, $clust) = split("\t", $line);
		my ($num, $ind) = $orignum =~ /^(\d+)\.(\d+)$/ if $orignum =~ /^\d+\.\d+$/;
		$ind = "" if not defined $ind;
		LOG($outLog, "\t\torignum=$orignum, num=$num, ind=$ind\n") if $linecount5 < 5;
		my $numlen = $end - $beg;
		push(@clust, $num) if not grep(/^$num$/, @clust);
		if (not defined $clustlen{$num} or (defined $clustlen{$num} and $clustlen{$num} < $numlen)) {
			$clust{$num} = "$num\.0\t$beg\t$end\t$y\t$y2\t$clust";
			$clustlen{$num} = $numlen;
		}
		LOG($outLog, date() . "Cannot parse num from line=$line\n") if not defined $num;
		my ($mid) = int(($end + $beg)/2+0.5);
		if (not defined $used{$num} or (defined $used{$num}{len} and $used{$num}{len} < ($end - $beg))) {
			$used{$num}{clust} = $clust;
			$used{$num}{len} = $end - $beg;
		}
		push(@{$cl{$clust}{beg}}, $beg);
		push(@{$cl{$clust}{mid}}, $mid);
		push(@{$cl{$clust}{end}}, $end);
	}
	close $in2;
	foreach my $num (keys %used) {
		my $clust = $used{$num}{clust};
		$cl{$clust}{total_read_unique}{$num} = 1;
	}

