#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G);
getopts("vd:n:G:");

#########
# BEGIN #
#########
#

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
my $uuid = getuuid();
my ($dist, $footPeakFolder) = ($opt_d, $opt_n);

# sanity check -n footPeakFolder
die "\nUsage: $YW$0$N $LGN-g gene$N $CY-n <footPeak's output folder (footPeak's -o)>$N\n\n" unless defined $opt_n and -d $opt_n;
my $outDir = "$footPeakFolder/FOOTCLUST/";
makedir($outDir);
die "Failed to create output directory $LCY$outDir$N!\n" unless -d $outDir;
makedir("$outDir/.TEMP");
die "Failed to create output directory $LCY$outDir/.TEMP/$N!\n" unless -d "$outDir/.TEMP";
makedir("$outDir/PNG/");
die "Failed to create output directory $LCY$outDir/PNG$N!\n" unless -d "$outDir/PNG/";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footClust.txt") or die "Failed to create outLog file $outDir/logFile_footClust.txt: $!\n";

# parse footLoop log file
($footPeakFolder) = getFullpath($footPeakFolder);
my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
my ($footLoopFolder);
my ($geneIndexFile);
open (my $infootPeak_logFile, "<", $footPeak_logFile) or die;
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
die "Cannot find footLoop Folder from logfile $footPeak_logFile\n" unless defined $footLoopFolder and -d $footLoopFolder;
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
LOG($outLog, "gene Index file = $geneIndexFile\n");

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
   my ($label2, $gene, $strand, $window, $thres, $type) = parseName($fileName1);# =~ /^(.+)_gene(.+)_(Unk|Pos|Neg)_(\d+)_(\d+\.?\d*)_(\w+)\.(PEAK|NOPK)$/;
   LOG($outLog, "Using label=$label2. Inconsistent label in filename $LCY$fileName1$N\nLabel from $footPeakFolder/.LABEL: $label\nBut from fileName: $label2\n\n") if $label ne $label2;
   $label = $label2;
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
		my ($num1, $num2, $num3) = $read =~ /^.*\.?m(\d+_\d+).+\/(\d+)\/(ccs|\d+_\d+)/;
		DIELOG($outLog, "\n\nERROR AT PARSING NUMBERS from read=$LGN$read$N\n\n") if not defined $num1 or not defined $num2 or not defined $num3;
		$num3 = "0" if $num3 eq "ccs";
		my $num = "$num1$num2$num3";
		$num =~ s/_//g;
		
#		my ($num) = $read =~ /^.+\/(\d+)\/ccs/; 
#		if (not defined $num) {
#			my ($num1, $num2) = $read =~ /^.+\/(\d+)\/(\d+\_\d+)/;
#			$num = "$num1\_$num2" if defined $num1 and defined $num2;
#		}
#		($num) = $read =~ /\.(\d+)$/ if not defined $num;
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
	.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\", \"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\", .libPaths()) )\nlibrary(labeling)\nlibrary(ggplot2)\nlibrary(reshape2)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(RColorBrewer)
	
	set.seed(420)
	setwd(\"$outDir\");
	library(ggplot2)
	df = read.table(\"$outDir/.TEMP/$fullName1.temp\",header=T,sep=\"\\t\",row.names=1,colClasses=c(\"factor\",\"integer\",\"integer\"))
	lenz = length(unique(paste(df\$beg,df\$end)))
	k.best = 0
	if (lenz > 20) {
		set.seed(420)
		for (i in 2:30) {
			k = kmeans(df,i)
			mylen = length(k[k\$withinss/1e6 > 3])
			if (mylen == 0) {
				k.best = i
				break
			}
		}
		if (k.best == 0) {
			lenz = 5
		} else { lenz = k.best}
		df\$cluster = k\$cluster
	} else if (lenz > 5) {
		k = kmeans(df,2)
		df\$cluster = k\$cluster
	} else {
		df\$cluster = 1
	}
	#dm = kmeans(df,lenz,nstart=20)
	df = df[order(df\$cluster, as.integer(df[,1]/200), as.integer(df[,2]/200)),]
	df\$y = seq(1,dim(df)[1])
	df\$ymax = seq(2,dim(df)[1]+1)
	colnames(df) = c(\"x\",\"xmax\",\"clust\",\"y\",\"ymax\")
	df2 = as.data.frame(aggregate(df[,c(1,2)],by=list(df\$clust),function(x)mean(x,trim=0.05)))
	colnames(df2) = c(\"clust\",\"x2\",\"y2\")
	df2 = df2[order(as.integer(df2\$x2/200), as.integer(df2\$y2/100)),]
	df2\$clust2 = seq(1,dim(df2)[1])
	df\$id = rownames(df)
	df = as.data.frame(merge(df,df2[,c(1,4)]))
	df\$clust = df\$clust2
	mylen=500;#as.integer((df\$xmax-df\$x)/300)*100
#	mylen[mylen < 100] = 100
#	mylen[mylen > 500] = 500
	df = df[order(df\$clust, as.integer(df\$x/mylen), as.integer(df\$xmax/mylen)),]
	df\$y = seq(1,dim(df)[1])
	df\$ymax = seq(2,dim(df)[1]+1)
	df[,c(1,6)] = df[,c(6,1)]
	colnames(df)[c(1,6)] = colnames(df)[c(6,1)]
	df = df[,-7]
	png(\"$outDir/PNG/$fullName1.clust.png\",height=(dim(df)[1])*5,width=max(df\$xmax)+10)
	ggplot(df, aes(x,y)) + geom_rect(aes(xmin=x,ymin=y,xmax=xmax,ymax=ymax,fill=as.factor(clust))) + theme_bw() + 
	theme(panel.grid=element_blank()) + coord_fixed(ratio=10)
	dev.off()
	write.table(df,\"$outDir/.TEMP/$fullName1.clustz\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
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
	my %clust; my %clustlen; my @clust; my $clustheader = "";
	open (my $in2, "<", "$outDir/.TEMP/$fullName1.clustz") or DIELOG($outLog, date() . __LINE__ . "\tFailed to read from $outDir/.TEMP/$fullName1.clustz: $!\n");
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
	open (my $outz2, ">", "$outDir/.TEMP/$fullName1.clust") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outDir/.TEMP/$fullName1.clust: $!\n");
	print $outz2 "$clustheader\n";
	foreach my $num (@clust[0..@clust-1]) {
		print $outz2 "$clust{$num}\n";
	}
	close $outz2;

	my ($ampliconLength) = $genes{uc($gene)};
	DIELOG($outLog, date() . " Failed to get amplicon length frmo gene=" . uc($gene) . ": $!\n") if not defined $ampliconLength;
	# process clust: get average beg/end point
	makedir("$outDir/CLUST_GENOME/") if not -d "$outDir/CLUST_GENOME";
	my $outgenome = "$outDir/CLUST_GENOME/$fullName1.clust.bed";
	$outgenome =~ s/.local/.genome/i;
	open (my $out2, ">", "$outDir/.TEMP/$fullName1.clust.bed") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outDir/.TEMP/$fullName1.clust.bed: $!\n");
	open (my $out3, ">", "$outgenome") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outgenome: $!\n");
	print $out2 "#gene\tbeg\tend\tcluster\ttotal_peak\tstrand\n";
	my ($coorCHR, $coorBEG, $coorEND, $coorSTRAND) = ($coor{uc($gene)}{chr}, $coor{uc($gene)}{beg}, $coor{uc($gene)}{end}, $coor{uc($gene)}{strand});
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
		my $total_read = scalar(@{$cl{$clust}{beg}});
		my $total_read_unique = scalar(keys %{$cl{$clust}{total_read_unique}});
		my $begGENOME = $coorBEG + $begOfALL;
		my $endGENOME = $coorBEG + $endOfALL;
		print "$coorCHR\t$begGENOME\t$endGENOME\t$gene.$clust\t$total_read_unique.$total_read\t$coorSTRAND\n";
		print $out3 "$coorCHR\t$begGENOME\t$endGENOME\t$gene.$clust\t$total_read_unique.$total_read\t$coorSTRAND\n";
		print $out2 "$gene\t$begOfALL\t$endOfALL\t$clust.WHOLE\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfBEG\t$endOfBEG\t$clust.BEG\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfMID\t$endOfMID\t$clust.MID\t$total_read_unique.$total_read\t$strand\n";
		print $out2 "$gene\t$begOfEND\t$endOfEND\t$clust.END\t$total_read_unique.$total_read\t$strand\n";
	}
	close $out2;
	LOG($outLog, "$YW ::: fastaFromBed -fi $faFile.footClustfastafile -bed $outDir/.TEMP/$fullName1.clust.bed -fo $outDir/.TEMP/$fullName1.clust.fa -s ::: $N\n","NA");
	system("fastaFromBed -fi $faFile.footClustfastafile -bed $outDir/.TEMP/$fullName1.clust.bed -fo $outDir/.TEMP/$fullName1.clust.fa -s > $outDir/.TEMP/$fullName1.clust.bed.LOG 2>&1") == 0 or LOG($outLog, date() . " Failed to fastaFromBed -fi$LCY $faFile.footClustfastafile$N -bed$LGN $outDir/.TEMP/$fullName1.clust.bed$N -fo$LPR $outDir/.TEMP/$fullName1.clust.fa$N -s >$YW $outDir/.TEMP/$fullName1.clust.bed.LOG$N 2>\&1: $!\n");

	LOG($outLog, date() . "$LCY\tAnalyzing Sequence Content from$N $outDir/.TEMP/$fullName1.clust.fa\n");
	open (my $faIn, "$outDir/.TEMP/$fullName1.clust.fa") or DIELOG($outLog, "Failed to read from $outDir/.TEMP/$fullName1.clust.fa: $!\n");
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
