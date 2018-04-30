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

my ($dist, $footPeakFolder) = ($opt_d, $opt_n);

# sanity check -n footPeakFolder
die "\nUsage: $YW$0$N $LGN-g gene$N $CY-n <footPeak's output folder (footPeak's -o)>$N\n\n" unless defined $opt_n and -d $opt_n;
($footPeakFolder) = getFullpath($footPeakFolder);
my $outDir = "$footPeakFolder/FOOTCLUST/";
makedir($outDir);
die "Failed to create output directory $LCY$outDir$N!\n" unless -d $outDir;
makedir("$outDir/.TEMP");
die "Failed to create output directory $LCY$outDir/.TEMP/$N!\n" unless -d "$outDir/.TEMP";
makedir("$outDir/PNG/");
die "Failed to create output directory $LCY$outDir/PNG$N!\n" unless -d "$outDir/PNG/";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footClust.txt") or die "Failed to create outLog file $outDir/logFile_footClust.txt: $!\n";

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
foreach my $input1 (sort @local_peak_files) {
	$input1_count ++;
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
	$strand = $strand eq "Neg" ? "-" : "+";
	# check peak file. Skip if there's less than 10 peaks
	my ($total_line) = `wc -l $input1` =~ /^(\d+)/;
	if ($total_line <= 10) {
		LOG($outLog, "\n--------------\n\n$LGN$input1_count$N. ${LRD}NEXTED $N: $input1 ${LCY}because total peaks is less than 10 $N($LGN$total_line$N)\n\n");
		$files_log .= "$label\t$gene\t$strand\t$window\t$thres\t$type\t${LRD}Skip$N\ttotal_peak_all=$total_line\ttotal_read_unique=-1\ttotal_peak_used=-1\tPEAKS_LOCAL=PEAKS_LOCAL/$fullName1\tPEAK_FILE=NA\n";
		next;
	}
	$files_log .= "$label\t$gene\t$strand\t$window\t$thres\t$type\t${LGN}Ran$N\ttotal_peak_all=$total_line";

	# Actual parse of peak file
	my ($linecount, $total_peak_used, $total_peak_all, %data) = (0,0,0);
	open (my $in1, "sort -k1,1 -k2,2n -k3,3n $input1|") or LOG($outLog, date() . "Cannot read from $input1: $!\n") and die;
	LOG($outLog, "\n--------------\n\n$LGN$input1_count$N. ${LCY}RUNNING$N: $input1\n");
	LOG($outLog, date() . "$LCY\tInfo$N: pcb=$label,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type\n");
	LOG($outLog, date() . "$LCY\tExample Lines$N:\n");
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my ($read, $beg, $end) = split("\t", $line);
		my ($num) = $read =~ /^.+\/(\d+)\/ccs/; 
			($num) = $read =~ /\.(\d+)$/ if not defined $num;
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
		push(@{$cl{$clust}{mid}}, $mid);
		push(@{$cl{$clust}{end}}, $end);
	}
	close $in2;
	
	# process clust: get average beg/end point
	open (my $out2, ">", "$outDir/.TEMP/$fullName1.clust.bed") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $outDir/.TEMP/$fullName1.clust.bed: $!\n");
	print $out2 "#gene\tbeg\tend\tcluster\ttotal_peak\tstrand\n";
	foreach my $clust (sort {$a <=> $b} keys %cl) {
		my $BEG   = int(tmm(@{$cl{$clust}{beg}})+0.5);
		my $BEGSD = int(tmmsd(@{$cl{$clust}{beg}})+0.5);
		my ($beg0, $end0) = ($BEG - $BEGSD, $BEG + $BEGSD);
		my $MID   = int(tmm(@{$cl{$clust}{mid}})+0.5);
		my $MIDSD = int(tmmsd(@{$cl{$clust}{mid}})+0.5);
		my ($beg1, $end1) = ($MID - $MIDSD, $MID + $MIDSD);
		my $END   = int(tmm(@{$cl{$clust}{end}})+0.5);
		my $ENDSD = int(tmmsd(@{$cl{$clust}{end}})+0.5);
		my ($beg2, $end2) = ($END - $ENDSD, $END + $ENDSD);
		my $total = @{$cl{$clust}{beg}};
		print $out2 "$gene\t$beg0\t$end2\t$clust.WHOLE\t$total\t$strand\n";
		print $out2 "$gene\t$beg0\t$end0\t$clust.BEG\t$total\t$strand\n";
		print $out2 "$gene\t$beg1\t$end1\t$clust.MID\t$total\t$strand\n";
		print $out2 "$gene\t$beg2\t$end2\t$clust.END\t$total\t$strand\n";
	}
	close $out2;
	open (my $faIn2, "<", $faFile) or die;
	open (my $faOut, ">", "$faFile.footClustfastafile") or die;
	my $fasta2 = new FAlite($faIn2);
	while (my $entry = $fasta2->nextEntry()) {
		my $def = uc($entry->def);
		my $seq = $entry->seq;
		print $faOut "$def\n$seq\n";
	}
	close $faIn2;
	close $faOut;
	system("/bin/rm $faFile.footClustfastafile.fai") if -e "$faFile.footClustfastafile.fai";
	LOG($outLog, "$YW ::: fastaFromBed -fi $faFile.footClustfastafile -bed $outDir/.TEMP/$fullName1.clust.bed -fo $outDir/.TEMP/$fullName1.clust.fa -s ::: $N\n");
	system("fastaFromBed -fi $faFile.footClustfastafile -bed $outDir/.TEMP/$fullName1.clust.bed -fo $outDir/.TEMP/$fullName1.clust.fa -s");

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
