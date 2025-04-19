#!/usr/bin/perl
	
use strict; use warnings FATAL=>'all'; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_t $opt_R $opt_D $opt_0 $opt_J $opt_F $opt_f);
getopts("vd:n:G:t:RD:0J:Ff");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib/';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/softwares/';
}

use myFootLib;
use FAlite;

my $homedir = $ENV{"HOME"};
my $usage = "\n$YW$0$N -n $LCY<footPeakFolder>$N\n\n";

my ($peakFile, $faFile, $label, $footLoopFolder, $geneIndexFile, $totalpeakFiles, $forcerun) = @ARGV;
my ($footPeakFolder, $max_parallel_run, $clustThreshold, $windowDiv, $dist) = ($opt_n, $opt_J, $opt_t, $opt_D, $opt_d);
die $usage  unless @ARGV >= 6 and defined $opt_n and -e $opt_n;

my @footPeakFolderShort = split("/", $footPeakFolder);
my ($footPeakFolderShort) = @footPeakFolderShort[@footPeakFolderShort-1];
my $toggleRstrand = defined $opt_R ? "Yes" : "No";

my $genewantprint = defined $opt_G ? $opt_G : "<N/A>";
my $forceprint = defined $opt_F ? $opt_F : "FALSE";
my @names = qw(peakFile faFile label footLoopFolder geneIndexFile totalpeakFiles);
print "\n$LGN############ footClust sbatch 2.pl #############$N\n";
print "\n$YW$0$N \\\n";
print "$LCY-n footPeakFolder$N \\ $LGN# -n genewantprint$N\n";
print "$LCY-J max_parallel_run$N \\ $LGN# -J max_parallel_run$N\n";
print "$LCY-t $clustThreshold$N \\ $LGN# -t clustThreshold$N\n";
print "$LCY-D $windowDiv$N \\ $LGN# -D windowDiv$N\n";
print "$LCY-d $dist$N \\ $LGN# -d dist$N\n";
print "$LCY-G $genewantprint$N \\ $LGN# -G genewant$N\n" if defined $opt_G;
print "$LCY-F$N \\ $LGN# force$N\n" if defined $opt_F;
for (my $i = 0; $i < @names; $i++) {
	print "$LCY$ARGV[$i]$N \\ $LGN#$names[$i]$N\n";
}
print "\n\n";

my $date = getDate();
my $uuid = getuuid();
my ($peakFolder, $peakFilename) = getFilename($peakFile, "folderfull");
my $outDir = getFullpath("$footPeakFolder/FOOTCLUST/");

# establish log file
my $footClust_logFile = "$footPeakFolder/.footClust_sbatch_2/$peakFilename\_footClust_logFile.txt";

open (my $outLog, ">", $footClust_logFile) or print "Failed to create outLog file $footClust_logFile: $!\n" and exit 1;
LOG($outLog, "$footClust_logFile\n");

my %coor;
open (my $inGeneIndexFile, "<", $geneIndexFile) or LOG($outLog, "Failed to read from $geneIndexFile: $!\n") and exit 1;
while (my $line = <$inGeneIndexFile>) {
   chomp($line);
   my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
   $coor{uc($gene)}{chr} = $chr;
   $coor{uc($gene)}{beg} = $beg;
   $coor{uc($gene)}{end} = $end;
   $coor{uc($gene)}{strand} = $strand;
}
close $inGeneIndexFile;

# Parse fasta file
my %genes;
open (my $faIn2, "<", $faFile) or die;
my $fasta2 = new FAlite($faIn2);
while (my $entry = $fasta2->nextEntry()) {
   my $def = uc($entry->def);
   my $seq = $entry->seq;
   my $length = length($seq);
   $def =~ s/>//;
   $genes{uc($def)} = $length;
}
close $faIn2;

foreach my $key1 (sort keys %genes) {
	LOG($outLog, date() . "${LPR}Check genehash$N: key=$LCY$key1$N val=$LGN$genes{$key1}$N\n");
	last;
}
LOG($outLog, "\n\n");

my @local_peak_files = ($peakFile);

####################
# Processing Input #
####################
my $files_log = "#FOLDER=$footPeakFolder\n";
my $local_peak_file_count = -1;
my $curr_gene = -1;
my $type_count = -1;
my @Rscript;

for (my $i = 0; $i < @local_peak_files; $i++) {
	my $local_peak_file = $local_peak_files[$i];
	if (defined $opt_G) {
		if ($local_peak_file !~ /$opt_G/i) {
			next;
		}
		else {
			LOG($outLog, date() . "${LPR}footClust_sbatch_2.pl$N: $LGN$i$N ${LGN}Processing$N $LCY$local_peak_file$N as it contain $LGN-G $opt_G$N\n");
		}
	}
	else {
		LOG($outLog, date() . "${LPR}footClust_sbatch_2.pl$N: $LGN$i$N ${LGN}Processing$N $LCY$local_peak_file$N\n");
	}

	# remove double // from folder
	$local_peak_file =~ s/[\/]+/\//g;
	my ($folder1, $fileName1) = getFilename($local_peak_file, "folderfull");
	my ($fullName1) = getFilename($local_peak_file, "full");

	# get gene and strand from file name
	my $parseName = parseName($fileName1);
   my ($label2, $gene, $strand, $window, $thres, $type) = @{$parseName->{array}};
   LOG($outLog, "\n------------- $YW WARNING!!! $N -------------\n\nUsing label=$label2. Inconsistent label in filename $LCY$fileName1$N\nLabel from $footPeakFolder/.LABEL: $label\nBut from fileName: $label2\n\n--------------- $YW WARNING!!! $N ---------------\n\n") if $label ne $label2;
   $label = $label2;
	my ($coorCHR, $coorBEG, $coorEND, $coorSTRAND) = ($coor{uc($gene)}{chr}, $coor{uc($gene)}{beg}, $coor{uc($gene)}{end}, $coor{uc($gene)}{strand});
	LOG($outLog, "gene=$LCY$gene$N\n");
	DIELOG($outLog, date() . "Fatal error at footClust.pl file=$local_peak_file, gene=$LCY$gene$N main gene index strand (\$coorSTRAND) is $coorSTRAND (has to be + or -)\n") if $coorSTRAND !~ /^[\+\-]$/;
	my $Rstrand = $coorSTRAND eq "+" ? "Pos" : "Neg";

	$thres *= 100 if $thres < 1;
	$gene   = uc($gene);

	if ($curr_gene ne $gene) {
		LOG($outLog, "\n$YW -------------- Processing $LGN$gene$N ------------- \n");
		$local_peak_file_count ++;
		$type_count = -1;
	}
	$type_count ++;
	$curr_gene = $gene;
	$strand = $strand eq "Neg" ? "-" : "+";

	# check number of peak
	my ($total_line) = `wc -l $local_peak_file` =~ /^(\d+)/;

	# Actual parse of peak file
	my ($linecount, $total_peak_used, $total_peak_all, %data) = (0,0,0);
	open (my $in1, "sort -k1,1 -k2,2n -k3,3n $local_peak_file|") or LOG($outLog, date() . "Cannot read from $local_peak_file: $!\n") and die;
	LOG($outLog, "\n");
	LOG($outLog, date() . "$LGN$local_peak_file_count.$type_count$N ${LCY}RUNNING$N: label=$LPR$label2$N gene=$LGN$gene$N strand=$LCY$strand$N window=$LGN$window$N thres=$LCY$thres$N type=$LGN$type$N input1=$local_peak_file\n");
	LOG($outLog, date() . "$LCY\tInfo$N: pcb=$label,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type\n");
	LOG($outLog, date() . "$LCY\tExample Lines$N:\n");
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my ($read, $beg, $end) = split("\t", $line);
		my ($id) = parse_readName($read, $outLog);
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
	die "Cannot get filenameshort of $local_peak_file ($LCY$fullName1$N)\n" unless defined $fileNameShort;
	my $peakFile = "$footPeakFolder/.CALL/$fileNameShort.PEAK.out";
	die "Cannot find PEAK file of $local_peak_file ($LCY$peakFile$N)\n" unless -e $peakFile;

	# Calculate average beg/mid/end and write a temp bed file for each peak
	LOG($outLog, date() . "$LCY\tInfo 2$N: total_peak_all=$LGN$total_peak_all$N,total_read_unique=$LGN$total_read_unique$N,total_peak_used=$LGN$total_peak_used$N\n");
	my $tempFile1_Peak = "$outDir/.TEMP/.$fullName1.clust.temp1";
	my $tempFile2_Clust = "$outDir/.TEMP/.$fullName1.clust.temp2";
	my $tempFile2_PNG = "$outDir/PNG/$fullName1.clust.png";

	if ((keys %data) == 0) {
		LOG($outLog, "footClust: No peak in $tempFile1_Peak so not writtn!\n");
		next;
	}
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
#	library(Cairo)
	set.seed(420)
	print(\"$tempFile1_Peak\")
	df.main = read.table(\"$tempFile1_Peak\",header=T,sep=\"\\t\",colClasses=c(\"factor\",\"integer\",\"integer\",\"integer\")) # changed df to df.main
	df.main.total_peak = unique(df.main[,4])

	last_cluster = 0

   dftemp = df.main
   dftemp\$len = dftemp\$xend - dftemp\$xbeg
   dftemp2 = aggregate(dftemp\$len, by=list(dftemp\$id), max)
   colnames(dftemp2) = c(\"id\",\"len\")
   dftemp = merge(dftemp,dftemp2,by=\"id\")
   dftemp = dftemp[dftemp\$len.x == dftemp\$len.y,]
   dftemp = dftemp[dftemp\$len.x == dftemp\$len.y,]
   df = subset(dftemp,select=c(-len.x,-len.y))
   df\$total_peak = 1

	df2 = df[,-1]
	lenz = length(unique(paste(df\$xbeg,df\$xend)))
	k.best = 0
	if (lenz >= 20) {
		set.seed(420)
		for (i in 2:10) {
			k = kmeans(df2,i)
			mylen = length(k[k\$withinss/1e6 > $clustThreshold])
			if (mylen == 0) {
				k.best = i
				break
			}
		}
		if (k.best == 0) {
			lenz = 5
			k.best = 3
		} else if (k.best >= 10) {
			k.best = 10
		} else {
			lenz = k.best
		}
		k = kmeans(df2,k.best)
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
	write.table(df,\"$tempFile2_Clust\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
	png(\"$tempFile2_PNG\",height=min(20000,(dim(df)[1])*5),width=max(df\$xend)+10)
	ggplot(df, aes(xbeg,ybeg)) + geom_rect(aes(xmin=xbeg,ymin=ybeg,xmax=xend,ymax=yend,fill=as.factor(clust))) + theme_bw() + 
	theme(panel.grid=element_blank()) + coord_fixed(ratio=10)
	dev.off()
	";


	open (my $outR, ">", "$tempFile1_Peak.R") or DIELOG($outLog, date() . __LINE__ . "Failed to write to $tempFile1_Peak.R: $!\n");
	print $outR $Rscript;
	close $outR;
	if (not -e "$tempFile2_Clust" or defined $opt_F or defined $opt_f) {
		system("Rscript $tempFile1_Peak.R") == 0 or DIELOG($outLog, date() . __LINE__ . "Failed to run Rscript $tempFile1_Peak.R: $!\n") and exit 1;
	}

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

	# process clust
	my %clust0;
	open (my $clustIn0, "<", $tempFile2_Clust) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $tempFile2_Clust: $!\n");
	while (my $line = <$clustIn0>) {
		chomp($line);
		if ($line =~ /xbeg\txend/) {next;}
		my ($id, $xbeg, $xend, $ybeg, $yend, $clust, $id2) = split("\t", $line);
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

	my ($ybeg, $yend) = (1,2);
	my $clust = 0; my $lastClustText = "NA";
	open (my $outClust0, ">", $tempFile2_Clust) or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $tempFile2_Clust: $!\n");

	if ($Rstrand eq "Neg" and $toggleRstrand eq "Yes") {
		foreach my $id (sort {$clust0{$b}{total_peak} <=> $clust0{$a}{total_peak} || $clust0{$b}{clust} <=> $clust0{$a}{clust} || $clust0{$b}{end} <=> $clust0{$a}{end} || $clust0{$b}{beg} <=> $clust0{$a}{beg} } keys %clust0) {
			my $clustText = join(",", @{$clust0{$id}{clustarr}});
			$clust ++ if $lastClustText ne $clustText;
			$lastClustText = $clustText;
		}
		$clust ++;
		$lastClustText = "NA";
		foreach my $id (sort {$clust0{$b}{total_peak} <=> $clust0{$a}{total_peak} || $clust0{$b}{clust} <=> $clust0{$a}{clust} || $clust0{$b}{end} <=> $clust0{$a}{end} || $clust0{$b}{beg} <=> $clust0{$a}{beg} } keys %clust0) {
			my $clustText = join(",", @{$clust0{$id}{clustarr}});
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
			my $clustText = join(",", @{$clust0{$id}{clustarr}});
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

	my $tempFile3_Bed = $tempFile2_Clust . ".bed";
	my $tempFile4_Fa  = $tempFile3_Bed . ".fa";
	open (my $out2, ">", "$tempFile3_Bed") or DIELOG($outLog, date() . __LINE__ . "\tFailed to write to $tempFile3_Bed: $!\n");
	open (my $outlocalBed, ">", $outlocalBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outlocalBedFile: $!\n");
	open (my $outgenomeBed, ">", $outgenomeBedFile) or DIELOG($outLog, date() . __LINE__ . "\tFailed to read form $outgenomeBedFile: $!\n");
	print $out2 "#gene\tbeg\tend\tcluster\ttotal_peak\tstrand\n";
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
	LOG($outLog, "$files_log\n");

	LOG($outLog, "\n\n");
	LOG($outLog, date() . "LOG FILES: \n");
	LOG($outLog, date() . "${LPR}tempFile1_Peak.R$N: $LCY$tempFile1_Peak.R$N\n");
	LOG($outLog, date() . "${LPR}tempFile1_Peak$N: $LCY$tempFile1_Peak$N\n");
	LOG($outLog, date() . "${LPR}tempFile2_Clust$N: $LCY$tempFile2_Clust$N\n");
	LOG($outLog, date() . "${LPR}tempFile2_Clust.fa$N: $LCY$tempFile2_Clust.fa$N\n");
	LOG($outLog, date() . "${LPR}tempFile2_PNG$N: $LCY$tempFile2_PNG$N\n");
	LOG($outLog, date() . "${LPR}tempFile3_Bed$N: $LCY$tempFile3_Bed$N\n");
	LOG($outLog, date() . "${LPR}tempFile3_Bed.LOG$N: $LCY$tempFile3_Bed.LOG$N\n");
	LOG($outLog, date() . "${LPR}tempFile4_Fa$N: $LCY$tempFile4_Fa$N\n");
	LOG($outLog, date() . "${LPR}LOG_FILESRAN$N: $LCY$outDir/.0_LOG_FILESRAN$N\n\n");
	LOG($outLog, date() . "${LPR}outlocalindivBedFile$N $LCY$outlocalindivBedFile$N\n");
	LOG($outLog, date() . "${LPR}outlocalindivClustFile$N $LCY$outlocalindivClustFile$N\n");
	LOG($outLog, date() . "${LPR}outgenomeindivBedFile$N $LCY$outgenomeindivBedFile$N\n");
	LOG($outLog, date() . "${LPR}outgenomeindivClustFile$N $LCY$outgenomeindivClustFile$N\n");
	LOG($outLog, date() . "${LPR}outlocalBedFile$N $LCY$outlocalBedFile$N\n");
	LOG($outLog, date() . "${LPR}outgenomeBedFile$N $LCY$outgenomeBedFile$N\n");
	LOG($outLog, "\n" . date() . "${LGN}SUCCESS!!!$N ${YW}footClust_sbatch_2.pl$N on $LCY$footPeakFolderShort$N\n\n");
}	
close $outLog;
