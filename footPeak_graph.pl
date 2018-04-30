#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_w $opt_g $opt_G $opt_v $opt_n); #v $opt_x $opt_R $opt_c $opt_t $opt_n);
getopts("n:vg:w:G:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";
my @version = `cd $footLoopDir && git log | head `;
my $version = "UNKNOWN";
foreach my $line (@version[0..@version-1]) {
   if ($line =~ /^\s+V\d+\.?\d*\w*\s*/) {
      ($version) = $line =~ /^\s+(V\d+\.?\d*\w*)\s*/;
   }
}
if (not defined $version or (defined $version and $version eq "UNKNOWN")) {
   ($version) = `cd $footLoopDir && git log | head -n 1`;
}
if (defined $opt_v) {
   print "\n\n$YW$0 $LGN$version$N\n\n";
   exit;
}

die "\nUsage: $YW$0$N [Optional: -G $LGN<gene to process>$N] -n$LCY <footPeak output directory>$N\n\n" if not defined $opt_n;
die "\nERROR: -n footPeak dir $LCY$opt_n$N doesn't exists!\n\nUsage: $YW$0$N -n <footPeak output directory>\n\n" if not -d $opt_n;

my $uuid = getuuid();
my $date = date();

main($opt_n);
sub main {
	my ($resDir) = @_;
	my %scp;
   makedir("$resDir/.CALL") if not -d "$resDir/\.CALL";
   makedir("$resDir/PEAKS_GENOME") if not -d "$resDir/PEAKS_GENOME";
   makedir("$resDir/PEAKS_LOCAL") if not -d "$resDir/PEAKS_LOCAL";
   makedir("$resDir/PNG") if not -d "$resDir/PNG";
   makedir("$resDir/PNG/PEAK") if not -d "$resDir/PNG/PEAK";
   makedir("$resDir/PNG/PEAKNEG") if not -d "$resDir/PNG/PEAKNEG";
   makedir("$resDir/PNG/NOPK") if not -d "$resDir/PNG/NOPK";
   makedir("$resDir/PNG/NOPKNEG") if not -d "$resDir/PNG/NOPKNEG";
   makedir("$resDir/PNG/ALL/") if not -d "$resDir/PNG/ALL/";
   makedir("$resDir/PDF") if not -d "$resDir/PDF";
   makedir("$resDir/PDF/PEAK") if not -d "$resDir/PDF/PEAK";
   makedir("$resDir/PDF/PEAKNEG") if not -d "$resDir/PDF/PEAKNEG";
   makedir("$resDir/PDF/NOPK") if not -d "$resDir/PDF/NOPK";
   makedir("$resDir/PDF/NOPKNEG") if not -d "$resDir/PDF/NOPKNEG";
   makedir("$resDir/PDF/ALL/") if not -d "$resDir/PDF/ALL/";
	my %files;
	my ($footPeak_logFile) = "$resDir/footPeak_logFile.txt";
	open (my $outLog, ">", "$resDir/footPeak_graph_logFile.txt") or die "\n\nFailed to write to $resDir/footPeak_graph_logFile.txt: $!\n\n";
	LOG($outLog, ">footPeak_graph.pl version $version\n");
	LOG($outLog, ">UUID: $uuid\n", "NA");
	LOG($outLog, ">Date: $date\n", "NA");
	LOG($outLog, ">Run script: $0 -n $opt_n\n", "NA");
	my %coor;
	my @lines = `cat $footPeak_logFile`;
	my ($label) = `cat $resDir/.LABEL`; chomp($label);
#	die "\n" if not -e "a";
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_logFile!\n\n") if not -e $footPeak_logFile;
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $resDir/.LABEL!\n\n") if not -e "$resDir/.LABEL";
	my ($thres, $window);
	foreach my $line (@lines) {
		chomp($line);
		if ($line =~ /^[ \t]+def=.+, coor=.+/) {
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
	   	if (defined $opt_G and $gene !~ /$opt_G/i) {
	   	   LOG($outLog, date() . " Skipped $LCY$gene$N as it doesn't contain $LGN-G $opt_G$N\n");
	   	   next;
	   	}
			$GENE = uc($GENE);
			DIELOG($outLog, "\n\ndied at processing $LCY$footPeak_logFile$N: can't parse index file def gene lqines\n\n$line\n\n") if not defined $STRAND;
			$coor{$GENE}{CHR} = $CHR;
			$coor{$GENE}{BEG} = $BEG;
			$coor{$GENE}{END} = $END;
			$coor{$GENE}{VAL} = $VAL;
			$coor{$GENE}{STRAND} = $STRAND eq "+" ? "Pos" : $STRAND eq "-" ? "Neg" : $STRAND =~ /^(Pos|Neg|Unk)$/ ? $STRAND : "Unk";
		}
		elsif ($line =~ /^-t thrshld\s+:/) {
			($thres) = $line =~ /^-t thrshld\s+:\s+(\-?\d+\.?\d*)$/;
		}
		elsif ($line =~ /^-w window\s+:/) {
			($window) = $line =~ /^-w window\s+:\s+(\-?\d+\.?\d*)$/;
		}
	}
	my %gene;
	my @types = qw(CH CG GH GC);
	foreach my $GENE (sort keys %coor) {
		my $mygene = $GENE;
		my $strand = $coor{$GENE}{STRAND};
		for (my $h = 0; $h < 4; $h++) {
			my $type = $types[$h];
			my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK";
			$files{$peakFile} = $mygene;
		}
	}
	my %Rscripts; 
	my $lastfile = -1; #debug
	LOG($outLog, "\n\nERROR: There is no file defined!\n") and die if (keys %files) == 0;#not defined $files_hash or (defined $files_hash and $files_hash !~ /HASH/);
	my $fileCount = 0;
	my $totalFile = (keys %files);
	my $lastGENE = -1;
	foreach my $file (sort keys %files) {
		$fileCount ++;
		my $GENE = $files{$file};
   	if (defined $opt_G and $file !~ /$opt_G/i) {
   	   LOG($outLog, date() . " Skipped $LCY$file$N as it doesn't contain $LGN-G $opt_G$N\n");
   	   next;
   	}
		LOG($outLog, "\n$YW -------- $fileCount/$totalFile Doing $GENE ---------$N\n\n") if $GENE ne $lastGENE;
		$lastGENE = $GENE;
		if (defined $opt_g and $GENE ne $opt_g) {
			LOG($outLog, date() . " $LCY Skipped $GENE$N (requested gene is $LGN$opt_g$N\n");
			next;
		} 
		if (defined $opt_w and $file !~ /$opt_w/) {
			LOG($outLog, date() . " Skipped $LCY$file$N (requested want is $LGN$opt_w$N\n");
			next;
		} 
#debug
		#next if $file !~ /CALM3.+Pos.+GC/;
#		last if $lastfile =~ /CALM3.+Pos.+/;
#debugend
		$lastfile = $file;
		my $STRAND = $coor{$GENE}{STRAND};
#		my $outPEAKS, ">", "$resDir/PEAKS_GENOME/$pk_filename.genome.bed")   or LOG($outLog, "\tFailed to write into $resDir/PEAKS_GENOME/$pk_filename.genome.bed: $!\n")  and exit 1;
#		open (my $outRPEAKS, ">", "$resDir/PEAKS_LOCAL/$pk_filename.local.bed") or LOG($outLog, "\tFailed to write into $resDir/PEAKS_LOCAL/$pk_filename.local.bed: $!\n") and exit 1;
		next if not defined $files{$file};
		my $peakFile = $file . ".out";
		my ($pk_filename) = getFilename($peakFile, 'full');
		$peakFile = "$resDir/.CALL/$pk_filename";
		my $nopkFile = $peakFile; $nopkFile =~ s/\.PEAK.out/.NOPK.out/;
		my $cluster_file = "$resDir/FOOTCLUST/.TEMP/$pk_filename";
		$cluster_file =~ s/.out$/.local.bed.clust/;
		my $kmer_file = "$resDir/FOOTCLUST/.TEMP/$pk_filename";
		$kmer_file =~ s/.out$/.local.bed.clust.kmer/;
	#	LOG($outLog, "$file");
	#	if (not -e $cluster_file) {
	#		LOG($outLog, "1=$LRD$cluster_file$N,");
	#	}
	#	if (-e $cluster_file) {
	#		LOG($outLog, "1=${LGN}$cluster_file$N,");
	#	}
	#	if (not -e $kmer_file) {
	#		LOG($outLog, "2=$LRD$kmer_file$N,");
	#	}
	#	if (-e $kmer_file) {
	#		LOG($outLog, "2=${LGN}$kmer_file$N,");
	#	}
	#	LOG($outLog, "\n");
#		my $nopkFile = $file . ".out"; 
		my $totpeak = -e $peakFile ? linecount($peakFile) : 0;
		my $totnopk = -e $nopkFile ? linecount($nopkFile) : 0;
		my $bedFile = "$resDir/PEAKS_LOCAL/$pk_filename.local.bed";
#		my $bedFile = $peakFile . ".RPEAKS"; 
		$bedFile =~ s/.out.local.bed/.local.bed/;
		my ($type) = $peakFile =~ /(CH|CG|GH|GC)/;
		LOG($outLog, date() . " $LGN$fileCount/$totalFile$N: Parsing into R script: gene=$LPR$GENE$N, strand=$GN$STRAND$N, type=$YW$type$N, peak=$LGN$totpeak$N, nopeak=$LRD$totnopk$N)\n");
		for (my $p = 0; $p < 2; $p ++) {
			my $currFile = $p == 0 ? $peakFile : $nopkFile;
			my $curr_cluster_file = $cluster_file;
			$curr_cluster_file =~ s/PEAK/NOPK/ if $p != 0;
			LOG($outLog, date() . " file #$YW$p.$N $LCY$currFile$N\n");
			LOG($outLog, "\t\tCurrfile           = $LCY$currFile$N
\t\tcurr_cluster_file  = $LCY$curr_cluster_file$N
\t\tkmer_File          = $LPR$kmer_file$N
\t\tbedFile            = $BU$bedFile$N
\t\ttotpeak = $LGN$totpeak$N, nopk = $LGN$totnopk$N
","NA");
			my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
			my ($label3, $gene3, $strand3, $window3, $thres3, $type3) = parseName($currFilename);
			#$currFilename =~ s/\.0_orig_//i;
			
			my $pngoutDir;
			if (($STRAND eq "+" or $STRAND eq "Pos") and $strand3 eq "Pos" and $type3 eq "CH" and $label3 =~ /PCB([1-9]$|10|12)/i) {
				$pngoutDir = $p == 0 ? "PEAK/" : "NOPK/";
			}
			elsif (($STRAND eq "-" or $STRAND eq "Neg") and $strand3 eq "Neg" and $type3 eq "GH" and $label3 =~ /PCB([1-9]$|10|12)/i) {
				$pngoutDir = $p == 0 ? "PEAK/" : "NOPK/";
			}
			elsif (($STRAND eq "+" or $STRAND eq "Pos") and $strand3 eq "Neg" and $type3 eq "GH" and $label3 =~ /PCB([1-9]$|10|12)/i) {
				$pngoutDir = $p == 0 ? "PEAKNEG/" : "NOPKNEG/";
			}
			elsif (($STRAND eq "-" or $STRAND eq "Neg") and $strand3 eq "Pos" and $type3 eq "CH" and $label3 =~ /PCB([1-9]$|10|12)/i) {
				$pngoutDir = $p == 0 ? "PEAKNEG/" : "NOPKNEG/";
			}
			elsif (($STRAND eq "+" or $STRAND eq "Pos") and $strand3 eq "Pos" and $type3 eq "CG" and $label3 =~ /PCB1[1345]/i) {
				$pngoutDir = $p == 0 ? "PEAK/" : "NOPK/";
			}
			elsif (($STRAND eq "-" or $STRAND eq "Neg") and $strand3 eq "Neg" and $type3 eq "GC" and $label3 =~ /PCB1[1345]/i) {
				$pngoutDir = $p == 0 ? "PEAK/" : "NOPK/";
			}
			elsif (($STRAND eq "+" or $STRAND eq "Pos") and $strand3 eq "Neg" and $type3 eq "GC" and $label3 =~ /PCB1[1345]/i) {
				$pngoutDir = $p == 0 ? "PEAKNEG/" : "NOPKNEG/";
			}
			elsif (($STRAND eq "-" or $STRAND eq "Neg") and $strand3 eq "Pos" and $type3 eq "CG" and $label3 =~ /PCB1[1345]/i) {
				$pngoutDir = $p == 0 ? "PEAKNEG/" : "NOPKNEG/";
			}
			else {
				$pngoutDir = "ALL/";
			}
			my $resDir2 = getFullpath($resDir);
			my $pngout = "$resDir/PNG/$pngoutDir$currFilename.png";
			my $pdfout = "$resDir/PDF/$pngoutDir$currFilename.pdf";
			my $lenpdfout = "$resDir/PDF/$pngoutDir$currFilename\_length.pdf";
			$scp{"scp mitochi\@crick.cse.ucdavis.edu:$resDir2/PNG/$pngoutDir$currFilename.png ./"} = 1;
			#LOG($outLog, "!HERE STRAND=$STRAND strand=$strand type=$type label=$label pngoutDir = $pngoutDir, p = $p, PNGOUT = $pngout\n");
#			next if not -e $currFile;
#			next if linecount($currFile) <= 1;
			my $Rscript = ".libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\", \"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\", .libPaths()) )\nlibrary(labeling)\nlibrary(ggplot2)\nlibrary(reshape2)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(RColorBrewer)\n";
			
			#if (-e $currFile and linecount($currFile) > 10) {
			my $totread = $totpeak + $totnopk;
			if (not -e $currFile or (-e $currFile and linecount($currFile) <= 10)) {
#				$Rscript .= "png(\"$pngout\",250,250)\nplot(seq(1,10),seq(1,10))\ntext(5,5,labels=c(\"PEAK = $totpeak, NOPK = $totnopk\"))\ndev.off()\n";
				$Rscript .= "png(\"$pngout\",1000,1000)\nplot(NA,xlim=c(1,100),ylim=c(1,100),xlab=NA,ylab=NA,bty=\"n\")\ntext(50,50,cex=3,labels=c(\"$currFilename\n\nPEAK = $totpeak / $totread\"))\ndev.off()\n";
			}
			else {
				my ($R) = Rscript($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout);

				# Read Table
				$Rscript .= $R->{readTable};

				# Cluster and Third Plot Cluster Graph
				if (-e $curr_cluster_file and -s $curr_cluster_file > 0) {
					$Rscript .= $R->{clusterFile};
				} 
				else {
					$Rscript .= $R->{noclusterFile};
				}

				# Read Bed File
				if (-e $bedFile and linecount($bedFile) > 0) {
					$Rscript .= $R->{peakbedFile};
				}

				# Main Plot
				$Rscript .= $R->{mainplot};

				# Main Plot Cluster Addition
				if (-e $curr_cluster_file and -s $curr_cluster_file > 0) {
					$Rscript .= $R->{mainplotClusterAddition};
				}
				# Main Plot Peak Bed Addition
				if (-e $bedFile and linecount($bedFile) > 0) {
					$Rscript .= $R->{mainplotPeakBedAddition};
				}
				
				# Second Plot & Conversion Graph
				$Rscript .= $R->{secondplotConversionGraph};

				# Add Third Plot and Do PNG
				$Rscript .= $R->{Scale};
#				$Rscript .= $R->{PDF};
				$Rscript .= $R->{PNG};
			}
			open (my $outRscript, ">", "$currFile.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.R: $!\n") and print $outLog $Rscript and next);
			print $outRscript $Rscript;
			$Rscripts{"$currFile.R"} = 1;
			close $outRscript;
		}
	}
	LOG($outLog, "\n\n$YW ----------------- Running R Scripts ------------------$N\n\n");
	$fileCount = 0;
	$totalFile = (keys %Rscripts);
	foreach my $outRscript (sort keys %Rscripts) {
		LOG($outLog, date() . " $LCY Skipped $outRscript$N (requested gene is $LGN$opt_g$N\n") and next if defined $opt_g and $outRscript !~ /$opt_g/;
		$fileCount ++;
		LOG($outLog, "\n" . date() . "$LGN$fileCount/$totalFile$N. Running $LCY$outRscript$N\n");
		LOG($outLog, date() . "\tR --vanilla --no-save < $outRscript > $outRscript.LOG 2>&1\n");
		my $RLOG = 0;
		$RLOG = system("R --vanilla --no-save < $outRscript > $outRscript.LOG 2>&1");
		my $prevRLOG = $RLOG;
		if ($RLOG ne 0) {
			if (not -e "$outRscript.LOG") {
				$RLOG = "\t$outRscript.LOG cannot be found!\n";
			}
			else {
				my @RLOG = `tail -n 5 $outRscript.LOG`;
				$RLOG = "\t" . join("\n\t", @RLOG);
			}
			LOG($outLog, date() . "\t${LRD}Failed$N to run_Rscript.pl $outRscript: $prevRLOG, LOG:\n$RLOG\n");
		}
		else {
			LOG($outLog, date() . "\t${LGN}Success$N on running run_Rscript.pl $LCY$outRscript$N\n");
		}
	}
#	LOG($outLog, date . "\tcd $resDir && run_Rscript.pl *MakeHeatmap.R\n");
#	system("cd $resDir && run_Rscript.pl *MakeHeatmap.R") if not defined $opt_x and defined $opt_R;
	LOG($outLog, "\n\n$YW ----------------- SCP PATHS ------------------$N\n\n");
	foreach my $file (sort keys %scp) {
		LOG($outLog, "$file\n");
	}
}

###############
# Subroutines #
###############


sub Rscript {
	my ($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout) = @_;
	my $R;
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	($bedFilename) =~ s/(CG|CH|GC|GH).+$/$1/;

# -------------------- $R->{readTable}
	$R->{readTable} = "	
# Read Table
df = read.table(\"$currFile\",skip=1,sep=\"\\t\")
colnames(df) = c(\"V1\",seq(1,dim(df)[2]-1))

# Sort table
if (dim(df)[1] < 1000) {
	h = hclust(dist(df[,-1]))
	df = df[h\$order,]
} else if (dim(df)[2] < 10) {
	mysum = apply(df[,-1],1,sum)
	df = df[order(mysum),]
} else {
	mybin = as.integer((dim(df)[2] - 1)/5)
	df.temp = df[,c(1,2)]
	currpos = 2
	for (i in 1:5) {
		df.temp[i+1] = apply(df[,seq(currpos, currpos+mybin)], 1, sum)
	}
	h = hclust(dist(df.temp[,-1]))
	df = df[h\$order,]
}

# Get read number ID
if (length(grep(\"ccs\",df\$V1)) > 0) {
	df\$id = gsub(\"^.+/([0-9]+)/ccs.\*\$\", \"\\\\1\", df\$V1, perl=T)
} else if (length(grep(\"\\\\.[0-9]+\$\",df\$V1)) > 0) {
	df\$id = gsub(\"^.+\\\\.([0-9]+)\$\", \"\\\\1\", df\$V1, perl=T)
} else {
	df\$id = df\$V1
}
						";

# -------------------- $R->{clusterFile}

	$R->{clusterFile} = "
# Cluster
clust = read.table(\"$curr_cluster_file\",header=T,sep=\"\\t\")
clust = clust[grep(\"^[0-9]+\\\\.[0-9]+\$\",clust\$id,perl=T,invert=T),]
clust\$y = seq(1,dim(clust)[1])
clust = subset(clust,select=c(\"id\",\"y\",\"clust\"))
clust2 = as.data.frame(aggregate(clust\$y,by=list(clust\$clust),min))
clust2\$max = aggregate(clust\$y,by=list(clust\$clust),max)\$x
colnames(clust2) = c(\"clust\",\"ymin\",\"ymax\")
clust2\$xmin = 1
clust2\$xmax = 30
clust2\$clust = clust2\$clust + 10
df3 = merge(df,clust,by=\"id\")
df3 = subset(df3,select=c(-y,-id,-V1))
df3clust = df3\$clust
df4 = as.data.frame(matrix(nrow=max(df3\$clust),ncol=dim(df3)[2]-1))
for (x in 1:max(df3\$clust)) {
	for (y in 1:(dim(df3)[2]-1)) {
		a = df3[df3\$clust == x,y]
		mymax = max(a)
		peak = length(a[a == 8 | a == 9])
		conv = length(a[a == 6 | a == 7])
		none = length(a[a == 4 | a == 5])
		if (mymax == 1) {
			df4[x,y] = 99
		} else if (peak != 0 & peak/2 >= conv) {#mymax == 8 | mymax == 9) {
			df4[x,y] = as.integer(length(a[a == 8 | a == 9]) / length(a) * 9+0.5)
		} else if (conv != 0 & conv > peak/2) {#mymax == 6 | mymax == 7) {
			df4[x,y] = as.integer(length(a[a == 6 | a == 7]) / length(a) * -9+0.5)
		} else if (mymax == 4 | mymax == 5) {
			df4[x,y] = 0
		}  else if (mymax == 0) {
			df4[x,y] = -99
		} else {
			df4[x,y] = 100
		}
	}
}

colnames(df4) = seq(1,dim(df4)[2])
df4\$clust = seq(1,max(df3clust))
df4\$y = df4\$clust
df3 = melt(df4,id.vars=c(\"y\",\"clust\"))
#df3 = melt(df3,id.vars=c(\"y\",\"clust\")); 
colnames(df3) = c(\"y\",\"clust\",\"x\",\"value\")
df3\$x = as.numeric(as.character(df3\$x))
df3\$clust = as.numeric(as.character(df3\$clust))
df3\$y = as.numeric(as.character(df3\$y))
df3\$value = as.numeric(as.character(df3\$value))
df5 = df3; df5[df5\$value > 10 | df5\$value < -10,]\$value = 0
df4 = data.frame(clust=seq(1,max(df5\$clust)),xmin=-1,xmax=-1,ymin=-1,ymax=-1)
for (i in (min(df5\$clust) : max(df5\$clust))) {
	if (length(df5[df5\$clust == i,]\$y) > 0) {
		df4\$xmin[i] = min(df5[df5\$value > 0 & df5\$clust == i,]\$x)-0.5
		df4\$xmax[i] = max(df5[df5\$value > 0 & df5\$clust == i,]\$x)+0.5
		df4\$ymin[i] = min(df5[df5\$value > 0 & df5\$clust == i,]\$y)-0.5
		df4\$ymax[i] = max(df5[df5\$value > 0 & df5\$clust == i,]\$y)+0.5
	} else {
		df4 = df4[-i,]
	}
}
df5 = data.frame(clust=seq(1,max(df5\$clust)),xmin=1,xmax=max(df5\$x),
ymin=seq(1,max(df5\$clust))-0.5,ymax=seq(1,max(df5\$clust))+0.5)
df3\$value = as.factor(df3\$value)
greens = rev(brewer.pal(9,\"Greens\"))
reds = brewer.pal(9,\"Reds\")
df3col=c(
\"-99\" = \"grey\",
\"-9\" = greens[1],
\"-8\" = greens[2],
\"-7\" = greens[3],
\"-6\" = greens[4],
\"-5\" = greens[5],
\"-4\" = greens[6],
\"-3\" = greens[7],
\"-2\" = greens[8],
\"-1\" = greens[9],
\"0\" = \"cornsilk\",
\"1\" = reds[1],
\"2\" = reds[2],
\"3\" = reds[3],
\"4\" = reds[4],
\"5\" = reds[5],
\"6\" = reds[6],
\"7\" = reds[7],
\"8\" = reds[8],
\"9\" = reds[9],
\"99\" = \"white\")
clust = subset(clust,select=c(\"id\",\"y\"))
df = merge(df,clust,by=\"id\")
df = subset(df,select=-id)
df3\$x = as.numeric(as.character(df3\$x))
df3\$y = as.numeric(as.character(df3\$y))
df3\$value  = as.numeric(as.character(df3\$value))
df4\$clust2 = df4\$clust + 10

p3 = ggplot(df3,aes(x,y)) +
	geom_tile(aes(fill=as.factor(value))) +
	geom_rect(data=df5,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x=xmin,y=ymin),
	size=0.8,fill=rgb(0,0,0,alpha=0),color=rgb(0,0,0,alpha=0.25)) +
	geom_rect(data=df4,aes(xmin=xmin,xmax=xmax,ymin=ymin,
		ymax=ymax,x=xmin,y=ymin,color=as.factor(df4\$clust2)),
		size=1,fill=rgb(0,0,0,alpha=0)) +
	geom_rect(data=df4,aes(fill=as.factor(df4\$clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax)) +
	geom_text(data=df4,aes(x=10,y=(ymin+ymax)/2,label=clust2-10),hjust=0,size=5) +
	theme_bw() + theme(legend.position=\"none\") + coord_cartesian(ylim=c(-1,6)) +
	scale_fill_manual(values=c(df3col,\"11\"=\"#e41a1c\",\"12\"=\"#377eb8\",
	\"13\"=\"#4daf4a\",\"14\"=\"#984ea3\",\"15\"=\"#ff7f00\")) +
	scale_color_manual(values=c(df3col,\"11\"=\"#e41a1c\",\"12\"=\"#377eb8\",
	\"13\"=\"#4daf4a\",\"14\"=\"#984ea3\",\"15\"=\"#ff7f00\")) +
	scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank())

";

# -------------------- $R->{noclusterFile}
	$R->{noclusterFile} = " 

# Cluster: NO CLUSTER
df\$y = seq(1,dim(df)[1])\ndf=subset(df,select=-id)

";

# -------------------- $R->{peakbedFile}

	$R->{peakbedFile} = "

# Peak BedFile
bed = read.table(\"$bedFile\",sep=\"\\t\")
bed = merge(subset(df,select=c(\"V1\",\"y\")),bed,by=\"V1\")

";

# -------------------- $R->{mainplot}
	$R->{mainplot} .= "
# Main Plot Part 1
dm = melt(df,id.vars=c(\"V1\",\"y\"))
dm\$variable = as.numeric(as.character(dm\$variable))

p =	ggplot(dm,aes(variable,y)) +  
		geom_tile(aes(fill=as.factor(value))) + 
	theme_bw() + theme(legend.position=\"none\") + 
	scale_fill_manual(
		values=c( 
			\"0\"=\"grey\",
			\"1\"=\"white\",
			\"4\"=\"cornsilk\",
			\"5\"=\"cornsilk\",
			\"6\"=\"green4\",
			\"7\"=\"seagreen4\",
			\"8\"=\"red4\",
			\"9\"=\"maroon4\",
			\"11\"=\"#e41a1c\",
			\"12\"=\"#377eb8\",
			\"13\"=\"#4daf4a\",
			\"14\"=\"#984ea3\",
			\"15\"=\"#ff7f00\"
		)
	) +
	scale_x_continuous(expand = c(0,0)) + 
	scale_y_continuous(expand = c(0,0)) +
	theme(
		line = element_blank(),
		axis.text = element_blank(),
		axis.title = element_blank()
	) + 
	ggtitle(paste(\"(peak=\",$totpeak,\"; nopk=\",$totnopk,\")\",sep=\"\"))

";


# -------------------- $R->{mainplotClusterAddition}
	$R->{mainplotClusterAddition} = "

# Main Plot Cluster Addition
p = 
	p + 
	geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
	geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=clust-10),hjust=0,size=15)

";

# -------------------- 	$R->{mainplotPeakBedAddition}
	$R->{mainplotPeakBedAddition} = "

# Main Plot Peak Bed Addition
if (length(bed) != 0 & dim(bed)[1] > 0) {
	bed\$variable=1
	p = 
		p + 
		geom_rect(data=bed,aes(xmin=V2,xmax=V3,ymin=y-0.5,ymax=y+0.5),
			size=0.5,fill=rgb(0,0,0,alpha=0),color=rgb(1,0,0,alpha=0.25))

	bed.q = quantile(bed\$V3-bed\$V2,probs=c(0.25,0.5,0.75))
	bed.maxy = max(bed\$V3-bed\$V2)
	if (bed.maxy <= 3000) {
		bed.maxy= 3000
	}
	p.bed = ggplot(bed,aes(x=1,y=V3-V2)) + 
			  annotate(geom=\"segment\",x=1,xend=1,y=0,yend=bed.maxy,lty=2,size=0.25,color=\"grey\") +
           annotate(geom=\"segment\",x=0.95,xend=1.05,y=bed.q[1],yend=bed.q[1],lty=2,size=0.5,color=\"black\") +
           annotate(geom=\"segment\",x=0.95,xend=1.05,y=bed.q[2],yend=bed.q[2],lty=1,size=1,color=\"black\") +
           annotate(geom=\"segment\",x=0.95,xend=1.05,y=bed.q[3],yend=bed.q[3],lty=2,size=0.5,color=\"black\") +
           annotate(geom=\"point\",x=1,y=bed.q[1],size=1,color=\"black\") +
           annotate(geom=\"point\",x=1,y=bed.q[2],size=2,color=\"black\") +
           annotate(geom=\"text\",x=1.1,y=bed.q[2]+50,label=as.character(bed.q[2]),size=4,color=\"black\") +
           annotate(geom=\"point\",x=1,y=bed.q[3],size=1,color=\"black\") +
			  geom_violin(fill=rgb(1,1,1,0)) +
			  theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\",
					axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
			  ylab(\"Peak length (bp)\") + xlab(\"$bedFilename\") +
			  coord_cartesian(ylim=c(-100,bed.maxy))

	pdf(\"$lenpdfout\",width=7,height=7)
	grid.arrange(p.bed)
	dev.off()
}

";


# -------------------- $R->{secondplotConversionGraph}
	my $peakminVal = $currFile =~ /PEAK.out$/ ? 8 : 6;
	$R->{secondplotConversionGraph} = "

# Calculate % Conversion
df2 = subset(df,select=c(-V1,-y));
df2[df2 < $peakminVal] = 0; df2[df2 >= 8] = 1
df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean))
if (dim(df2[df2\$y > 0,])[1] > 15) {
	df2 = df2[df2\$y > 0,]
	df2\$x = as.numeric(as.character(df2\$x));
	df2\$y = as.numeric(as.character(df2\$y));
	df2\$x2 = df2\$x
	df2\$y2 = df2\$y
	for (i in 1:(dim(df2)[1]-10)) {
		a = df2[df2\$x >= df2[i,]\$x & df2\$x <= df2[i+10-1,]\$x,]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2[i,]\$y2 = mean(a\$y)
			df2[i,]\$x2 = mean(a\$x)
		}
	}
	mins = seq(1,as.integer(df2[1,]\$x2)-1,10)
	maxs = seq(max(df2\$x2),dim(df)[2]-2,10)
	df2 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2)
	df2 = rbind(df2,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2 = data.frame(x=seq(1,dim(df)[2]), y=0, x2=seq(1,dim(df)[2]), y2=0);
}

# P2 % Conversion XY Plot
p2 = 
	ggplot(df2,aes(x2,y2)) + geom_point(aes(x=x,y=y),size=1) + geom_line(color=rgb(1,0,0,alpha=1)) + theme_bw()+
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{Scale} = "

ratio1 = dim(df)[1] / (dim(df)[1] + 31.25 + 26.5625)
ratio2 = 31.25      / (dim(df)[1] + 31.25 + 26.5625)
ratio3 = 26.5625    / (dim(df)[1] + 31.25 + 26.5625)
mynrow = 2
totalheight = dim(df)[1] + 31.25
totalwidth  = dim(df)[2] * 0.125
totalratio  = c(ratio1, ratio2)
	if (file.exists(\"$curr_cluster_file\") & file.info(\"$curr_cluster_file\")\$size > 0) {
		totalheight = totalheight + 26.5625
		totalratio = c(totalratio, ratio3)
		mynrow = 3
	}

# Scaling
myscale = 4
totalheight = totalheight * myscale
totalwidth = totalwidth * myscale
totalratio = totalratio * myscale

";

	$R->{PNG} = "

# PNG
png(\"$pngout\",width=totalwidth,height=totalheight)
if (mynrow == 3) {
	grid.arrange(p,p2,p3,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p,p2,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

";

	$R->{PDF} = "

# PDF
pdf(\"$pdfout\",width=totalwidth,height=totalheight)
if (mynrow == 3) {
	grid.arrange(p,p2,p3,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p,p2,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

";

	return $R;
}

sub make_total_hash {
	my $total;
	my @types = qw(CH CG GH GC);
	foreach my $type (@types) {
		$total->{$type}{peak} = 0; 
		$total->{$type}{nopk} = 0;
		$total->{$type}{total} = 0;
	}
	return($total);
}
sub make_heatmap {
	

}





sub parse_peak {
	my ($ARG, $bad, $minDis, $minLen, $outLog) = @_;
	my ($name, $isPeak, $mygene, $type, $strand, @val) = split("\t", $ARG);
	my %bad = %{$bad} if defined $bad;
	my $name_want = "AIRN_PFC66_FORWARD.16024";#CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/16024/ccs";
	shift(@val) if $val[0] eq "";
#	for (my $i = 0; $i < @val; $i++) {
#		if ($val[$i] !~ /^[456789]$/) {print "."} else {print "$val[$i]";}
#		LOG($outLog, date() . "\n") if $i != 0 and ($i+1) % 100 == 0;
#	}
#	LOG($outLog, date() . "\n");
#	0001234000
#	0123456789
#	len=10, e1=len(e1), e2=10-len(e2)
	my $peaks;
	my %peak; $peak{curr} = 0; #my $edge = 0; my $edge2 = 0; my $zero = 0; my $edge1 = 0;
	my $Length = @val; 
	my $print = "name=$name, isPeak = $isPeak, Total length = $Length\n";
	my ($edge1) = join("", @val) =~ /^(0+)[\.1-9A-Za-z]/;
	$edge1 = defined $edge1 ? length($edge1) : 0;
	my ($edge2) = join("", @val) =~ /[\.1-9A-Za-z](0+)$/;
	$edge2 = defined $edge2 ? @val-length($edge2) : @val;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		if ($i % 100 == 0) {$print .= "\n$YW" . $i . "$N:\t";}
		if ($val[$i] =~ /[89]/) {
			$peak{beg} = $i if $peak{curr} == 0;
			$print .= "${LPR}$val[$i]$N" if $peak{curr} == 0;
			$print .= "${LRD}$val[$i]$N" if $peak{curr} == 1;
			$peak{curr} = 1;
		}
		elsif ($val[$i] =~ /[23]/) {
			$peak{end} = $i+1;
			push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
			undef $peak{beg}; undef $peak{end};
			$peak{curr} = 0;
			$val[$i] =~ tr/23/89/;
			$print .= "${LPR}$val$N";
		}
		else {
			$print .= "EDGE1" if $i == $edge1;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[57]$/;
			$print .= "." if $val[$i] eq 1;
			$print .= "x" if $val[$i] eq 0;# and $i < $edge1;
			$print .= "EDGE2" if $i == $edge2 - 1;
		}
	}
	my (%nopk, @peak);
	$strand = $type =~ /^C/ ? 0 : $type =~ /^G/ ? 16 : 255;
	my %peak2;
	$print .= "\n";
	print "$print" if $name_want eq $name;
#	LOG($outLog, date() . "\nDoing $YW$name$N\n" if $name eq $name_want;#"SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
	if (defined $peak{peak}) {
		foreach my $peak (sort @{$peak{peak}}) {
			my ($beg, $end) = split("-", $peak);
#			LOG($outLog, date() . "$name: $beg to $end\n" if $name eq 77011 or $name eq "$name_want");
			my $checkBad = 0;
#			if (not defined $bad->{$strand}) {
#				push(@peak, "$beg-$end");
#				push(@{$peak2{peak}}, $peak);
#			}
#			next if not defined $bad->{$strand};
			foreach my $begBad (sort keys %{$bad->{$strand}}) {
				my $endBad = $bad->{$strand}{$begBad};
#			foreach my $begBad (sort keys %bad) {
#				my $endBad = $bad->{$strand}{$begBad};
#				LOG($outLog, date() . "\t$beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "$name_want");
				if ($beg >= $begBad and $beg <= $endBad and $end >= $begBad and $end <= $endBad) {
#					for (my $m = $beg; $m <= $end; $m++) {
#						$nopk{$m} = 1;
#					}
					LOG($outLog, date() . "\t\t$LGN YES$N peak=$beg-$end, bad=$begBad-$endBad\n") if $isPeak eq "PEAK";# if $name eq "$name_want");
					$checkBad = 1; last;
				}
			}
			if ($checkBad != 1) {
				next if not defined $bad->{$strand};
				foreach my $begBad (sort keys %{$bad->{$strand}}) {
					my $endBad = $bad->{$strand}{$begBad};
#					LOG($outLog, date() . "$name: $beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "$name_want");
					#if (($beg >= $begBad and $beg <= $endBad) or ($end >= $begBad and $end <= $endBad)) {
						my @valz = @val;
						my ($goodC, $badC) = (0,0);
						for (my $m = $beg; $m < $end; $m++) {
							if ($m >= $begBad and $m <= $endBad) {
								$badC ++ if $valz[$m] =~ /[2389]/;
							}
							else {
								$goodC ++ if $valz[$m] =~ /[2389]/;
							}
						}
						if ($goodC < 5 and $badC >= 9) {
#							LOG($outLog, date() . "\t$LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
#							LOG($outLog, date() . "\t$YW$name$N $LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n");
							$checkBad = 1; last;
						}
#						else {
#							#LOG($outLog, date() . "\t$LGN OKAY!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
#						}
#					LOG($outLog, date() . "\t\t$LGN YES$N\n" if $name eq "$name_want");
					#}
				}
			}
#			LOG($outLog, date() . "\t$name checkbad = $checkBad\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
			
			if ($checkBad == 1) {
				$print .= "\tCheckBad; Peak Not: $LRD$peak$N\n";
#				LOG($outLog, date() . "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak bad : $LRD$peak$N\n" if $name eq "$name_want");
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
			elsif ($end - $beg < $minLen) {
				$print .= "\tend-$end < $minLen; Peak Not: $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
			elsif ($beg < $edge2 - 100 and $end > 100 + $edge1) {
#				LOG($outLog, date() . "something wrong\n";# if $name eq "$name_want");
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Used: $LGN$peak$N\n";
				push(@peak, "$beg-$end");
				push(@{$peak2{peak}}, $peak);
			}
			else {
				$print .= "\tbeg=$beg, end=4end, peak Not: $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopk{$j} = 1;
				}
			}
		}
	}
	my $totalpeak = scalar(@peak);
	my @val2 = @val;
	if ($totalpeak > 0) {
		for (my $i = 0; $i < @val; $i++) {
			my $val = $val[$i];
			$val2[$i] = $val;
			if ($val =~ /^(8|9)$/ and defined $nopk{$i}) { # Peak Converted CpG or CH
				$val2[$i] = 7 if $val eq 9;
				$val2[$i] = 6 if $val eq 8;
			}
		}
	}
	#die $print if $totalpeak > 1;
	$print .= "$name\t$totalpeak\n" if $isPeak eq "PEAK";
#	LOG($outLog, date() . "$print\n" if $isPeak eq "PEAK";# and $print =~ /; Peak Not/;# if $totalpeak == 1;# or $name eq "SEQ_100022") and exit 1;
#	exit 0 if $isPeak eq "PEAK";
	@val = @val2;
	$print .= "\n\nVAL2: Total length = $Length\n";
	($edge1) = join("", @val) =~ /^(0+)[\.1-9A-Za-z]/;
	$edge1 = defined $edge1 ? length($edge1) : 0;
	($edge2) = join("", @val) =~ /[\.1-9A-Za-z](0+)$/;
	$edge2 = defined $edge2 ? @val-length($edge2) : @val;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		if ($i % 100 == 0) {$print .= "\n$YW" . $i . "$N:\t";}
		if ($val[$i] =~ /[89]/) {
			$peak{beg} = $i if $peak{curr} == 0;
			$print .= "${LPR}$val[$i]$N" if $peak{curr} == 0;
			$print .= "${LRD}$val[$i]$N" if $peak{curr} == 1;
			$peak{curr} = 1;
		}
		elsif ($val[$i] =~ /[23]/) {
			$peak{end} = $i+1;
			push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
			undef $peak{beg}; undef $peak{end};
			$peak{curr} = 0;
			$val[$i] =~ tr/23/89/;
			$print .= "${LPR}$val$N";
		}
		else {
			$print .= "EDGE1" if $i == $edge1;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[57]$/;
			$print .= "." if $val[$i] eq 1;
			$print .= "x" if $val[$i] eq 0;# and $i < $edge1;
			$print .= "EDGE2" if $i == $edge2 - 1;
		}
	}
	$print .= "\n";
	print "$print" if $name_want eq $name;

	return ($name, \@val2, $totalpeak, $peak2{peak});
}

sub find_lots_of_C {
	my ($seqFile, $mygene, $outLog) = @_;#, $geneIndex, $box) = @_; #$geneIndexesFa;
	my %seq;
	LOG($outLog, date . "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n");
	open(my $SEQIN, "<", $seqFile) or LOG($outLog, date() . "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!") and exit 1;
	my $fasta = new FAlite($SEQIN);
	my %lotsOfC;

	while (my $entry = $fasta->nextEntry()) {
	   my $gene = uc($entry->def); $gene =~ s/^>//;
	   my $seqz = uc($entry->seq);
		next if $gene ne $mygene;
	   LOG($outLog, date . "\t\tgenez=$gene ($gene)\n");
	   
	
		my $minlen = 6;
	   my $seqz2 = $seqz;#join("", @{$seq{$gene}{seq}});
	   while ($seqz2 =~ /(C){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_C) = length($curr);
	      my ($next_C) = $next =~ /^(C+)[AGTN]*$/;
	      $next_C = defined $next_C ? length($next_C) : 0;
	      my ($beg_C) = defined $prev ? length($prev) : 0;
	      my ($end_C) = $curr_C + $next_C + $beg_C;
	      my $length = $curr_C + $next_C;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_C to $end_C ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n");
	      $lotsOfC{$gene} .= "C;$beg_C;$end_C,";
	   }
		$seqz2 = "";
	   $seqz2 =$seqz;# join("", @{$seq{$gene}{seq}});
	   while ($seqz2 =~ /(G){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_G) = length($curr);
	      my ($next_G) = $next =~ /^(G+)[ACTN]*$/;
	      $next_G = defined $next_G ? length($next_G) : 0;
	      my ($beg_G) = defined $prev ? length($prev) : 0;
	      my ($end_G) = $curr_G + $next_G + $beg_G;
	      my $length = $curr_G + $next_G;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_G to $end_G ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n");
	      $lotsOfC{$gene} .= "G;$beg_G;$end_G,";
	   }
		my $bad;
		if (defined $lotsOfC{$gene}) {
			$lotsOfC{$gene} =~ s/,$//;
			my @lotsOfC = split(",", $lotsOfC{$gene});
			foreach my $coor (@lotsOfC) {
				my ($nuc, $beg, $end) = split(";", $coor);
				my $strands = $nuc eq "C" ? 0 : $nuc eq "G" ? 16 : 255;
				$bad->{$strands}{$beg} = $end-1;
				LOG($outLog, date() . "$strands: coor=$coor, nuc=$nuc, beg=$beg, end=$end, beg-end=$beg-$bad->{$strands}{$beg}\n");
			}
			return $bad;
		}
	}
	return;
}
#foreach my $gene (keys %lotsOfC) {
#   $gene = uc($gene);
#   $lotsOfC{$gene} =~ s/;$//;
#   LOG($outLog, date() . "$gene\t$lotsOfC{$gene}\n");
#   my $beg2 = $geneIndex{$gene};
#   foreach my $lines (@{$box->{$gene}}) {
#      LOG($outLog, date() . "GENEZ = $gene, lines = $lines\n");
#   }
#   LOG($outLog, date() . "genez=$gene,beg=$beg2\n");
#}
#push


1;
__END__
# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

   # For nucleotide
# 10 = Nucleotide A
# 11 = Nucleotide C
# 12 = Nucleotide T
# 13 = Nucleotide G






__END__
END1
#my (@inputs, $input1);
#if ($input1 =~ /.orig$/) {
#}
#else {
#	(@inputs) = <$folders/*Pos50.orig>;
#	LOG($outLog, date() . "Must have 1 input only! (" . scalar(@inputs) . "):\n" . join("\n-", @inputs) . "\n") and exit 1 if @inputs != 1 and not defined $opt_x;
#	$input1 = defined $opt_x ? $folders : $inputs[0];
#	LOG($outLog, date() . "INPUT1=$input1\n");
#}

	
			else {
            $Rscript .= " 
		png(\"$pngout\",width=dim(df)[2]*2,height=dim(df)[1]*16 + 500)
grid.arrange(p,p2,ncol=1,nrow=2,heights=c(ratio1,ratio2));\ndev.off()\n";
			}


__END__
=comment
sub main {
	# From footPeak.pl: 
	# (($peakFilez, $seqFile, $gene, $minDis, $resDir, $minLen, $SEQ));
	my ($input1, $faFile, $mygene, $minDis, $resDir, $minLen, $SEQ) = @_;

	my @foldershort = split("\/", $resDir);
	my $foldershort = pop(@foldershort);
	print date() . "\nusage: $YW$0$N [-c to use cpg] $CY<CALM3_Pos_20_0.65_CG.PEAK>$N $CY<location with lots of C>$N\n\n" and exit 1 unless @_ == 7;
	print date() . "Input cannot be directry!\n" and exit 1 if -d $input1;
	($input1) = getFullpath($input1);
	my ($folder, $fileName) = getFilename($input1, "folderfull");

	makedir("$resDir/.CALL") if not -d "$resDir/\.CALL";
	makedir("$resDir/PEAKS_GENOME") if not -d "$resDir/PEAKS_GENOME";
	makedir("$resDir/PEAKS_LOCAL") if not -d "$resDir/PEAKS_LOCAL";
	makedir("$resDir/PNG") if not -d "$resDir/PNG";
	makedir("$resDir/PNG/PEAK") if not -d "$resDir/PNG/PEAK";
	makedir("$resDir/PNG/PEAKNEG") if not -d "$resDir/PNG/PEAKNEG";
	makedir("$resDir/PNG/NOPK") if not -d "$resDir/PNG/NOPK";
	makedir("$resDir/PNG/NOPKNEG") if not -d "$resDir/PNG/NOPKNEG";
	makedir("$resDir/PNG/ALL/") if not -d "$resDir/PNG/ALL/";
	makedir("$resDir/PDF") if not -d "$resDir/PDF";
	makedir("$resDir/PDF/PEAK") if not -d "$resDir/PDF/PEAK";
	makedir("$resDir/PDF/PEAKNEG") if not -d "$resDir/PDF/PEAKNEG";
	makedir("$resDir/PDF/NOPK") if not -d "$resDir/PDF/NOPK";
	makedir("$resDir/PDF/NOPKNEG") if not -d "$resDir/PDF/NOPKNEG";
	makedir("$resDir/PDF/ALL/") if not -d "$resDir/PDF/ALL/";
	open (my $outLog, ">>", "$resDir/footLoop_addition_logFile.txt") or die;

	
	my @coor = split("\t", $SEQ->{$mygene}{coor});
	my (%pk, %Rscripts, %files);
	my $total = make_total_hash();

	my $label = "";
	if (-e "$resDir/.LABEL") {
	   ($label) = `cat $resDir/.LABEL`;
	   chomp($label);
	}
	else {
		DIELOG($outLog, "Failed to parse label from .LABEL in $resDir/.LABEL\n");
	}
	my ($label2, $gene, $strand, $window, $thres, $type) = parseName($fileName);# =~ /^(.+)_gene(.+)_(Unk|Pos|Neg)_(\d+)_(\d+\.?\d*)_(\w+)\.(PEAK|NOPK)$/;
	my $isPeak = $fileName =~ /\.PEAK/ ? "PEAK" : "NOPK";
	LOG($outLog, "Using label=$label2. Inconsistent label in filename $LCY$fileName$N\nLabel from $resDir/.LABEL: $label\nBut from fileName: $label2\n\n") if $label ne $label2;
	$label = $label2;

	if (defined $mygene) {
		LOG($outLog, date() . "footPeak.pl filename=$LCY$fileName$N, mygene=$mygene, input genez=$gene are not the same!\n") and return -1 if uc($mygene) ne uc($gene);;
	} else {$mygene = $gene;}
	
	my $bad = find_lots_of_C($faFile, $mygene, $outLog) if defined $faFile;

	LOG($outLog, date() . "$input1; Undefined mygene=$mygene=, strand=$strand=, window=$window=, thres=$thres=, type=$type=, isPeak=$isPeak=\n") and exit 1 if not defined $isPeak or not defined $window;
	LOG($outLog, date . "\n\nFolder $YW$folder$N: Processing files related to $LCY$input1$N\n");
	my @types = qw(CH CG GH GC);
	for (my $h = 0; $h < 4; $h++) {
		my $type = $types[$h];
		my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK";
		my $nopkFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.NOPK";
		$files{$peakFile} = 1;
		LOG($outLog, date . "h=$LGN$h\t$YW$peakFile\t$LCY$nopkFile\n$N");
	
		my ($folder1, $peakfileName) = getFilename($peakFile, "folderfull");
		my ($folder2, $nopkfileName) = getFilename($nopkFile, "folderfull");

		my $data;
		my ($linecount, $totalpeak, $totalnopk, $totalline) = (0,0,0,0);
		if (-e $nopkFile) {
			($totalline) = `wc -l $nopkFile` =~ /^(\d+)/;
			$linecount = 0;
			open (my $in1, "<", $nopkFile) or LOG($outLog, date() . "Cannot read from $nopkFile: $!\n") and exit 1;
			LOG($outLog, date . "\tProcessing NOPK file $LPR$nopkFile$N ($LGN$totalline$N lines)\n");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . "\tDone $totalnopk / $totalline\n") if $totalnopk % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, $bad, $minDis, $minLen, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data->{peak}}, $val) if $totalPeak > 0;
				push(@{$data->{nopk}}, $val) if $totalPeak == 0;
				$totalnopk ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in1;
		}
		my $peakCount = defined $data->{peak} ? @{$data->{peak}} : 0;
		my $nopkCount = defined $data->{nopk} ? @{$data->{nopk}} : 0;
		my $nopkPrint ="$folder2\t$nopkfileName\t$peakCount\t$nopkCount\t$totalnopk\t$totalline";
		LOG($outLog, date() . "$nopkPrint\n");
		$total->{$type}{peak}  += $peakCount;
		$total->{$type}{nopk}  += $nopkCount;
		$total->{$type}{total} += $totalnopk;
		
		if (-e $peakFile) {
			($totalline) = `wc -l $peakFile` =~ /^(\d+)/;
			$linecount = 0;
			open (my $in1, "<", $peakFile) or LOG($outLog, date() . "Cannot read from $peakFile: $!\n") and exit 1;
			LOG($outLog, date . "\tProcessing PEAK file $LPR$peakFile$N ($LGN$totalline$N lines)\n");
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1; #header
				LOG($outLog, date . "\tDone $totalpeak / $totalline\n") if $totalpeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, $bad, $minDis, $minLen, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data->{peak}}, $val) if $totalPeak > 0;
				push(@{$data->{nopk}}, $val) if $totalPeak == 0;
				$totalpeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in1;
		}
		$peakCount = defined $data->{peak} ? @{$data->{peak}} - $peakCount : 0;
		$nopkCount = defined $data->{nopk} ? @{$data->{nopk}} - $nopkCount : 0;
		my $peakPrint ="$folder1\t$peakfileName\t$peakCount\t$nopkCount\t$totalpeak\t$totalline";
		LOG($outLog, date() . "$peakPrint\n");
		$total->{$type}{peak}  += $peakCount;
		$total->{$type}{nopk}  += $nopkCount;
		$total->{$type}{total} += $totalpeak;

		if (defined $data->{peak}) {
			die if @{$data->{peak}} != $total->{$type}{peak};
			print "HERE: $folder1/$peakfileName.out\n";
			open (my $out1, ">", "$resDir/.CALL/$peakfileName.out") or LOG($outLog, date() . "Cannot write to $peakfileName.out: $!\n") and exit 1;
			foreach my $val (sort @{$data->{peak}}) {
				print $out1 "$val\n";
			}
			close $out1;
		}
		if (defined $data->{nopk}) {
			die if @{$data->{nopk}} != $total->{$type}{nopk};
			print "HERE: $folder1/$nopkfileName.out\n";
			open (my $out1, ">", "$resDir/.CALL/$nopkfileName.out") or LOG($outLog, date() . "Cannot write to $nopkfileName.out: $!\n") and exit 1;
			foreach my $val (sort @{$data->{nopk}}) {
				print $out1 "$val\n";
			}
			close $out1;
		}		
		LOG($outLog, date . "#Folder\tFile\tPeak\tnopk\tTotalRead\tTotalLineInFile\n") if $h == 0;
		LOG($outLog, date . "$peakPrint\n$nopkPrint\n");
	}
	my ($chr0, $beg0, $end0, $name0, $val0, $strand0) = @coor;
	my $STRAND = $strand0;
	die "Undefined beg or end at coor=\n" . join("\n", @coor) . "\n" if not defined $beg0 or not defined $end0;
	
#	system("footPeak_HMM.pl -n $resDir");

	foreach my $file (sort keys %pk) {
		my ($pk_filename) = getFilename($file, 'full');
		open (my $outPEAKS, ">", "$resDir/PEAKS_GENOME/$pk_filename.genome.bed")   or LOG($outLog, "\tFailed to write into $resDir/PEAKS_GENOME/$pk_filename.genome.bed: $!\n")  and exit 1;
		open (my $outRPEAKS, ">", "$resDir/PEAKS_LOCAL/$pk_filename.local.bed") or LOG($outLog, "\tFailed to write into $resDir/PEAKS_LOCAL/$pk_filename.local.bed: $!\n") and exit 1;
		my $currtype = $file =~ /_CH/ ? "CH" : $file =~ /_CG/ ? "CG" : $file =~ /_GH/ ? "GH" : $file =~ /_GC/ ? "GC" : "UNK";
		LOG($outLog, "\tFailed to determine type (CH/CG/GH/GC) of file=$file.PEAKS: $!\n") if $currtype eq "UNK";
		foreach my $name (sort keys %{$pk{$file}}) {
			foreach my $peak (sort @{$pk{$file}{$name}}) {
				my ($beg, $end) = split("-", $peak);
				my $end1 = $end + $beg0;
				my $beg1 = $beg + $beg0;
				my ($junk, $readname) = $name =~ /^(\w+)\.(.+)$/;
				print $outPEAKS  "$chr0\t$beg1\t$end1\t$name\t0\t$strand0\t$file\n";
				print $outRPEAKS "$name\t$beg\t$end\n";
			}
		}
		close $outPEAKS;
		close $outRPEAKS;
	}
	
	open (my $outLGENE, ">", "$resDir/.0_RESULTS\_$label\_gene$mygene\_$strand\_$window\_$thres.TXT");
	for (my $h = 0; $h < 4; $h++) {
		my $type = $types[$h];
		my $totalPeak = $total->{$type}{peak};
		my $totalNopk = $total->{$type}{nopk};
		$total->{$type}{peak} = $total->{$type}{total} == 0 ? 0 : int(1000 * $total->{$type}{peak} / $total->{$type}{total}+0.5)/10;
		$total->{$type}{nopk} = $total->{$type}{total} == 0 ? 0 : int(1000 * $total->{$type}{nopk} / $total->{$type}{total}+0.5)/10;
		my @folder = split("/", $resDir);
		my $foldershort = $folder[@folder-1];
		   $foldershort = $folder[@folder-2] if not defined ($foldershort) or (defined $foldershort and $foldershort =~ /^[\s]*$/);
		my $peakFile    = "$mygene\_$strand\_$window\_$thres\_$type.PEAK";
		print $outLGENE "#folder\tpeakFile\tGene\tStrand\ttotal\tpeak.perc\n" if $type eq "CH";
		print $outLGENE "$foldershort\t$peakFile\t$mygene\t$type\t$total->{$type}{total}\t$totalPeak\t$total->{$type}{peak}\n";
	}
	close $outLGENE;
	system("cat $resDir/.0_RESULTS\_$label\_gene$mygene\_$strand\_$window\_$thres.TXT");
	#my ($sampleName) = $folder =~ /\/?\d+_(m\d+_\d+)_\d+_\w+/;
	my ($sampleName) = $folder =~ /^.+(PCB[\_\-]*\d+)/i;
	if ($folder =~ /debarcode/) {
		my ($temp) = $folder =~ /_ccs_(\w+)/;
		$sampleName = $sampleName . "_$temp";
	}
	if (not defined $sampleName) {
		$sampleName = $foldershort;
	}
	system("footClust.pl -n $resDir -g $gene") == 0 or LOG($outLog, "Failed to run footClust.pl : $!\n");
	system("footClust2.pl -n $resDir -g $gene") == 0 or LOG($outLog, "Failed to run footClust2.pl : $!\n");
=cut

