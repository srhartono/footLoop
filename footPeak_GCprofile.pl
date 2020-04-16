#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_n $opt_i $opt_S $opt_G $opt_w $opt_F $opt_D);
getopts("vn:i:SG:w:FD");

#########
# BEGIN #
######### 

BEGIN {
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
	print "- Pushed $libPath into perl lib path INC\n\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
}

use myFootLib; use FAlite;

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

my $DEBUG = "NA" if not defined $opt_D;
my @treats = qw(gccont gcwskew purineskew atwskew);

################
# ARGV Parsing #
###############

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N -i $LGN<geneIndexFile.feature>$N -n $CY<footPeak output folder>$N

${LGN}Options:$N
-S: Skip the first 2 steps (e.g. if you've ran it before and would like to reuse the same files, as these usually don't change)
-G <gene>: Only run files with this gene in the name

";

die $usage unless defined $opt_n and defined $opt_i and -e $opt_i and -d $opt_n;

my ($indexFile, $footPeakFolder) = ($opt_i, $opt_n);
my ($genewant) = $opt_G if defined $opt_G;

# sanity check -n footPeakFolder

#($footPeakFolder) = getFullpath($footPeakFolder);
$footPeakFolder =~ s/\/$//;
my $footClustFolder = "$footPeakFolder/FOOTCLUST";
my $footKmerFolder  = "$footPeakFolder/KMER";
my $uuid = getuuid();
my ($user) = $homedir =~ /home\/(\w+)/;
my $date = date();


############
# LOG FILE #
############
open (my $outLog, ">", "$footPeakFolder/footPeak_GCprofile_logFile.txt") or die date() . ": Failed to create outLog file $footPeakFolder/footPeak_GCprofile_logFile.txt: $!\n";
LOG($outLog, ">$0 version $version\n");
LOG($outLog, ">UUID: $uuid\n", $DEBUG);
LOG($outLog, ">Date: $date\n", $DEBUG);
LOG($outLog, ">Run script: $0 -i $opt_i -n $opt_n\n", $DEBUG);


##########
# OUTDIR #
##########
my $resDir = "$footPeakFolder/GCPROFILE";
my ($resDirFullpath) = getFullpath($resDir);

makedir($resDir) if not -d $resDir;
makedir("$resDir/.TEMP") if not -d "$resDir/.TEMP";
makedir("$resDir/.TSV") if not -d "$resDir/.TSV";

my $OUTDIRS;
($OUTDIRS->{PDF}) = makeOutDir($resDirFullpath . "/PDF/");
foreach my $OUTDIR (keys %{$OUTDIRS->{PDF}}) {
	my $DIR = "$resDir/PDF/$OUTDIR/";
	next if not -d $DIR;
	my @pdf = <$DIR/*.pdf>;
	system("rm $DIR/*.pdf") if @pdf != 0;
	@pdf = <$DIR/ALL/*.pdf>;
	system("rm $DIR/ALL/*.pdf") if @pdf != 0;
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
	$gene{$gene}{strand} = $strand eq "+" ? "Pos" : $strand eq "-" ? "Neg" : $strand eq "." ? "Pos" : $strand;
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
my @bedFiles = <$footPeakFolder/FOOTCLUST/CLUST_GENOME/*.genome.bed.indiv.clust.bed>;
my @files;
my %coreFile;
my $bedFileCount = 0;
my %data;
LOG($outLog, "\n\n" . date() . "1. Getting cluster and preprocessing bed files! Folder:\n$LPR$footPeakFolder/FOOTCLUST/.TEMP$N\n\n");
foreach my $bedFile (sort @bedFiles) {
	if (defined $opt_G) {next if $bedFile !~ /$opt_G/};
	my ($bedFileName) = getFilename($bedFile, "full");
	my ($coreOutFile) = $bedFileName =~ /^(.+_CG|.+_CH|.+_GC|.+_GH).*$/;
	# process clustFile to get number of unique and non unique read with peak
	my ($clustFile) = $bedFile =~ /^(.+).bed$/;$clustFile =~ s/.indiv.clust/.clust/;
	open (my $clustFileIn, "<", $clustFile) or DIELOG($outLog, date() . "Failed to read from $clustFile: $!\n");
	while (my $line = <$clustFileIn>) {
		my ($chr, $beg, $end, $geneclust, $value, $strand) = split("\t", $line);
		my ($gene, $clust) = split(/\./, $geneclust);
		DIELOG($outLog, date() . "CAn't parse gene=$gene clust=$clust from geneclust=$geneclust line\n$line\nfile=$clustFile\n\n") if not defined $gene or not defined $clust;
		my ($totalreadunique, $totalread) = split(/\./, $value);
		$data{totalreadclust}{$coreOutFile}{$clust} = $totalread;
		$data{totalreaduniqueclust}{$coreOutFile}{$clust} = $totalreadunique;
		$data{totalread}{$coreOutFile} += $totalread;
		$data{totalreadunique}{$coreOutFile} += $totalreadunique;
#		die "$coreOutFile\n$clust\n$totalread\n$totalreadunique\n";
	}
	my @alphabets = qw(A B C W D E F);
	for (my $i = 0; $i < @alphabets; $i++) {
		for (my $j = 0; $j < @treats; $j++) {
			my $tsvFile = "$footPeakFolder/GCPROFILE/.TSV/$bedFileName\_100\_$alphabets[$i].temp.fa.$treats[$j].tsv";
			$coreFile{$tsvFile} = $coreOutFile;
			@files = (@files, $tsvFile);
		}
	}
	$bedFileCount ++;
	$cluster = get_cluster($bedFile, $coreOutFile, $cluster, $bedFileCount, $outLog);
	preprocess_bed($bedFile, $coreOutFile, $cluster, $outLog) if not defined $opt_S;
}

#########################
# Calculate Seq profile #
#########################

LOG($outLog, "\n\n" . date() . "2. Calculating Sequence profile (might take a couple minutes)\n\n");
if (not defined ($opt_S)) {
	my $res1 = ""; my $res2 = "";
	my @bedFiles = <$resDir/.TEMP/*_[ABCDEFW].temp>;
	my $fileCount = 0;
	foreach my $bedFile (@bedFiles[0..@bedFiles-1]) {
		$fileCount ++;
		my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
		my $cmd = "fastaFromBed -fi $genomeFile -bed $bedFile -fo $bedFolder/$bedFilename.fa -s -name";
		LOG($outLog, date() . "$YW$fileCount.$N Creating seq Profile of $LCY$bedFile$N\n");
		LOG($outLog, date() . "\t- $cmd\n","NA");
		$res1 .= `$cmd`;

		$cmd = "$footLoopScriptsFolder/lib/counter_cpg_indiv.pl -w 200 -s 1 -o $resDir/.TSV/ -A $bedFolder/$bedFilename.fa";
		LOG($outLog, date() . "\t- $cmd\n","NA");
		$res2 .= `$cmd`;
	}
	open (my $outRes1, ">", "$resDir/.TEMP/fastaFromBed.LOG") or LOG($outLog, "Cannot write to $resDir/.TEMP/fastaFromBed.LOG: $!\n");
	print $outRes1 "$res1\n";
	close $outRes1;
	open (my $outRes2, ">", "$resDir/.TEMP/counter_cpg_indiv.LOG") or LOG($outLog, "Cannot write to $resDir/.TEMP/counter_cpg_indiv.LOG: $!\n");
	print $outRes2 "$res2\n";
	close $outRes2;
}

######################
# Processing tsvFile #
######################


LOG($outLog, "\n\n" . date() . "3. Processing " . scalar(@files) . " files in $LCY$resDir$N\n");
my @header = ("wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$resDir/$N!\n" and exit if (@files == 0);
foreach my $tsvFile (sort @files) {
	LOG($outLog, "- Doing $LGN$tsvFile$N\n");
	my ($coreOutFile) = $coreFile{$tsvFile};#$tsvFile =~ /^(.+)_100_.\.temp.fa.\w+.tsv$/;

	$tsvFile =~ s/\/+/\//g;
	my ($folder1, $fileName1) = getFilename($tsvFile, "folderfull");
	if (defined $opt_G) {
		next if $fileName1 !~ /$opt_G/;
	}

   # get gene and strand from file name
   my $parseName = parseName($fileName1);
   my ($label2, $gene, $strand, $window, $thres, $type) = @{$parseName->{array}};

	next if $fileName1 !~ /^PCB/;
	my @arr = $fileName1 =~ /PEAK.genome.bed.indiv.clust.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	for (my $i = 0; $i < @arr; $i++) {
		DIELOG($outLog, date() . "fileName=$fileName1 Undefined i=$i header=$header[$i] arr[i] undef\n") if not defined $arr[$i];
	}
	my ($WINDOW, $SAMPLE, $CALCTYPE) = @arr;
	my ($outName) = $fileName1 =~ /^(.+).PEAK.genome.bed/;
	$outName .= ".$CALCTYPE";
	$data{data}{$outName} = $parseName;
	$data{data}{$outName}{sample} = $SAMPLE;
	$data{data}{$outName}{wind2} = $WINDOW;
	$data{data}{$outName}{calctype} = $CALCTYPE;
	LOG($outLog, date() . "  - gene=$LGN$gene$N pos=$LGN$SAMPLE$N calctype=$LPR$CALCTYPE$N input=$outName$N\n",$DEBUG) if $CALCTYPE eq "skew";
	open (my $in1, "<", $tsvFile) or DIELOG($outLog, date() . " Cannot read from $tsvFile: $!\n");
	my $linecount = 0;
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		$linecount ++;
		my ($read, $value) = split("\t", $line);
		my ($gene, $id, $cluster) = split(/\./, $read);
		#my ($read_number) = $read =~ /\.(m\d+.+)$/;
		#my ($id) = parse_readName($read_number);
		#my $cluster = $c$cluster->{$coreOutFile}{$id}{clust}; $cluster = -1 if not defined $cluster;
		#die "line=$line\ngene=$gene, id=$id, clust=$cluster\n";
		$data{cluster}{$outName}{$read} = $cluster;
		$data{id}{$outName}{$read} = $id;
		$data{read}{$outName}{$read}{$SAMPLE} = $value;
		$data{input}{$outName}{$SAMPLE} = $outName;
	}	
	close $in1;
}
open (my $out1, ">", "$resDir/RESULT.TSV") or DIELOG($outLog, date() . "Cannot write to $resDir/RESULT.TSV: $!\n");
foreach my $outName (sort keys %{$data{input}}) {
	my $WINDOW = $data{data}{$outName}{wind2};
	my $CALCTYPE = $data{data}{$outName}{calctype};
	my $SAMPLE = $data{data}{$outName}{sample};
	print $out1 "file\tid\tcluster\tread\twindow\ttype\tfeature";
	foreach my $sample (sort keys %{$data{input}{$outName}}) {
		print $out1 "\t$sample";
	}
	print $out1 "\toutfile\tflag\ttotalreadclust\ttotalreaduniqueclust\ttotalread\ttotalreadunique\n";
	last;
}

foreach my $outName (sort keys %{$data{read}}) {
	my ($coreOutFile) = $outName =~ /^(.+)\.\w+$/;
	my $WINDOW = $data{data}{$outName}{wind2};
	my $CALCTYPE = $data{data}{$outName}{calctype};
	my $SAMPLE = $data{data}{$outName}{sample};
	my $GENE = $data{data}{$outName}{gene};
	my $feature = $gene{$GENE}{feature}; 

   # get gene and strand from file name
   my $parseName = parseName($outName);
 	my $readStrand = $parseName->{strand};
 	my $rconvType  = $parseName->{type};
	my $geneStrand = $gene{$GENE}{strand};
	$feature = "FEATURE_UNKNOWN" if not defined $feature;
 	my $flag = getFlag($outName, $geneStrand, $readStrand, $rconvType);
	my $sampleoutDir = $resDirFullpath . "/PDF/$flag/";
	#print "$outName: gene=$GENE feature=$feature readStrand=$readStrand geneStrand=$geneStrand rconvtype=$rconvType $LCY$sampleoutDir$N\n";
	foreach my $read (sort keys %{$data{read}{$outName}}) {
		my $curr_cluster = $data{cluster}{$outName}{$read};
		my $totalreadclust = $data{totalreadclust}{$coreOutFile}{$curr_cluster};
		my $totalreaduniqueclust = $data{totalreaduniqueclust}{$coreOutFile}{$curr_cluster};
		my $totalread = $data{totalread}{$coreOutFile};
		my $totalreadunique = $data{totalreadunique}{$coreOutFile};
		DIELOG($outLog, date() . "Failed to get totalread from read=$read\n$outName\n$curr_cluster\n") if not defined $totalread or not defined $totalreadunique;
		print $out1 "$data{input}{$outName}{$SAMPLE}\t$data{id}{$outName}{$read}\t$data{cluster}{$outName}{$read}\t$read\t$WINDOW\t$CALCTYPE\t$feature";
		foreach my $sample (sort keys %{$data{read}{$outName}{$read}}) {
			print $out1 "\t$data{read}{$outName}{$read}{$sample}";
		}
		print $out1 "\t$sampleoutDir/$outName.GCPROFILE.$flag\t$flag\t$totalreadclust\t$totalreaduniqueclust\t$totalread\t$totalreadunique\n";
	}
}

GCprofile_Rscript($resDir, $outLog);

sub GCprofile_Rscript {
	my ($resDir, $outLog) = @_;
my $resDirFullpath = getFullpath($resDir);

my $RESULT = $resDir . "/RESULT.TSV";
my $LABEL = `cat $resDir/../.LABEL`; DIELOG($outLog, "Cannot find $resDir/../.LABEL!\n") if not defined $LABEL; chomp($LABEL);
my $LABEL2 = $LABEL;
$LABEL = $resDirFullpath . "/$LABEL";

LOG($outLog, "

If R dies for any reason, make sure you have these required R libraries:
- RColorBrewer v1.1-2
- gridExtra v2.3
- labeling v0.3
- reshape2 v1.4.3
- ggplot2 v3.1.0
- GMD v0.3.3

");


my $Rscript = "

library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)    

RESULT=\"$RESULT\";

dfmain = read.table(RESULT,header=T,sep=\"\\t\")
a = as.vector(sapply(dfmain,class))
a[2] = \"character\"
dfmain = read.table(RESULT,header=T,sep=\"\\t\",colClasses=a)
outfiles = unique(dfmain\$outfile)

# pdf all
mytypes = unique(as.character(dfmain\$type))
for (j in 1:length(mytypes)) {
	df = dfmain[dfmain\$type == mytypes[j],]
	dm = melt(df,id.vars=c(\"file\",\"id\",\"cluster\",\"read\",\"window\",\"type\",\"feature\",\"outfile\",\"flag\",\"totalreadclust\",\"totalreaduniqueclust\",\"totalread\",\"totalreadunique\"));
	dm\$variable = factor(dm\$variable,levels=c(\"A\",\"B\",\"C\",\"W\",\"D\",\"E\",\"F\"))
	dm\$gene = paste(dm\$file,dm\$feature)
	genes = unique(dm\$gene)
	dm\$group = paste(dm\$cluster)
	ylimsMin = -1
	ylimsMax = 1
	ylines = 0
	if(length(grep(\"dens\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylimsMin=0;
		ylimsMax=1.2;
		ylines=0.6;
	} else if (length(grep(\"cont\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylimsMin=0;
		ylimsMax=1;
		ylines=0.5;
	} else if (length(grep(\"wskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylimsMin=0.5;
		ylimsMax=-0.5;
		ylines=0;
	}
	
	if (length(grep(\"atwskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"AT Weighted Skew\"
	} else if (length(grep(\"atskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"AT Skew\"
	} else if (length(grep(\"gcskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"GC Skew\"
	} else if (length(grep(\"gcwskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"GC Weighted Skew\"
	} else if (length(grep(\"gccont\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"GC Content\"
	} else if (length(grep(\"cpgdens\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"CpG Density\"
	} else if (length(grep(\"purineskew\",mytypes[j],perl=T,ignore.case=T)) != 0) {
		ylabs = \"Purine Skew (G|A vs. C|T)\"
	} else {
		ylabs=mytypes[j]
	}
	myclust = unique(df\$cluster)
	uniquefeature = as.character(unique(df\$feature))
	uniqueflag = as.character(unique(df\$flag))

	for (k in 1:length(uniquefeature)) {
		for (l in 1:length(uniqueflag)) {
			mytitle   = paste(\"$LABEL2\.\",mytypes[j],\"\\n\",uniquefeature[k],\" \",uniqueflag[l],sep=\"\")
			outpdfall = paste(\"$LABEL\.\",mytypes[j],\".\",uniquefeature[k],\".\",uniqueflag[l],\".pdf\",sep=\"\")
			print(paste(j,\"Doing pdf\",outpdfall))
			temp = dm[dm\$feature == uniquefeature[k] & dm\$flag == uniqueflag[l],]
			if (length(temp) > 0) {
	
			pdf(outpdfall,width=7,height=7)
			mytotalreadunique = length(unique(temp\$id))
			mytotalread = length(unique(temp\$read))
			
			mytitle2 = paste(mytitle,\" (\",mytotalreadunique,\" unique reads, \",mytotalread,\" total peaks)\",sep=\"\")
			mytitle2 = paste(mytitle2,\"\\n\",ylabs,sep=\"\")

			p = ggplot(temp,aes(variable,value)) +
				geom_boxplot(aes(fill=variable),outlier.shape=NA) +
				theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin,ylimsMax)) +
				annotate(geom=\"segment\",x=0,xend=8,y=ylines,yend=ylines,lty=2) +
				ylab(ylabs) + xlab(\"Samples\") + ggtitle(mytitle2)

			print(p)
			dev.off()
			} else {print(\"\\tDoes not exist!\")}
		}
	}
}
for (j in 1:length(outfiles)) {
print(paste(j,\".\",sep=\"\"))
outfile=outfiles[j]
df = dfmain[dfmain\$outfile == outfiles[j],]
dm = melt(df,id.vars=c(\"file\",\"id\",\"cluster\",\"read\",\"window\",\"type\",\"feature\",\"outfile\",\"flag\",\"totalreadclust\",\"totalreaduniqueclust\",\"totalread\",\"totalreadunique\"));
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
types = as.character(df\$type[1])

ylimsMin = -1
ylimsMax = 1
ylines = 0
if(length(grep(\"dens\",types[1],perl=T,ignore.case=T)) != 0) {
	ylimsMin=0;
	ylimsMax=1.2;
	ylines=0.6;
} else if (length(grep(\"cont\",types[1],perl=T,ignore.case=T)) != 0) {
	ylimsMin=0;
	ylimsMax=1;
	ylines=0.5;
} else if (length(grep(\"wskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylimsMin=0.5;
	ylimsMax=-0.5;
	ylines=0;
}

if (length(grep(\"atwskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"AT Weighted Skew\"
} else if (length(grep(\"atskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"AT Skew\"
} else if (length(grep(\"gcskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"GC Skew\"
} else if (length(grep(\"gcwskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"GC Weighted Skew\"
} else if (length(grep(\"gccont\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"GC Content\"
} else if (length(grep(\"cpgdens\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"CpG Density\"
} else if (length(grep(\"purineskew\",types[1],perl=T,ignore.case=T)) != 0) {
	ylabs = \"Purine Skew (G|A vs. C|T)\"
} else {
	ylabs=types[1]
}
myclust = unique(df\$cluster)
uniquefile = as.character(unique(df\$file))
uniqueflag = as.character(unique(df\$flag))
mytitle = paste(uniquefile,\"\\n\",uniqueflag,sep=\"\")
for (i in 1:length(myclust)) {
	print(paste(j,\", Doing pdf of cluster i=\",i))
   temp = dm[dm\$cluster == myclust[i],]
	tempcount = length(unique(temp\$id))
	outfile2 = paste(outfile,\".cluster\",myclust[i],\".pdf\",sep=\"\")
	mytotalreadclust = as.character(unique(temp\$totalreadclust))
	mytotalreaduniqueclust = as.character(unique(temp\$totalreaduniqueclust))
	print(paste(i,tempcount,mytotalreadclust,mytotalreaduniqueclust))
	mytitle2 = paste(mytitle,\" cluster \",myclust[i],\" (\",mytotalreaduniqueclust,\" unique reads, \",mytotalreadclust,\" total peaks)\",sep=\"\")
	mytitle2 = paste(mytitle2,\"\\n\",ylabs,sep=\"\")
	pdf(outfile2,width=7,height=7)
   p = ggplot(temp,aes(variable,value)) +
      geom_boxplot(aes(fill=variable),outlier.shape=NA) +
      theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin,ylimsMax)) +
      annotate(geom=\"segment\",x=0,xend=8,y=ylines,yend=ylines,lty=2) +
      ylab(ylabs) + xlab(\"Samples\") + ggtitle(mytitle2)
   print(p)
   dev.off()
}

print(paste(j,\"Doing pdf\"))
outfile=paste(outfile,\".pdf\",sep=\"\")
	pdf(outfile,width=7,height=7)
   temp = dm;
	tempcount = length(unique(temp\$id))
	mytotalread = as.character(unique(temp\$totalread))
	mytotalreadunique = as.character(unique(temp\$totalreadunique))
	print(paste(\"Total = \",tempcount,mytotalread,mytotalreadunique))
	mytitle2 = paste(mytitle,\" (\",mytotalreadunique,\" unique reads, \",mytotalread,\" total peaks)\",sep=\"\")
	mytitle2 = paste(mytitle2,\"\\n\",ylabs,sep=\"\")
   p = ggplot(temp,aes(variable,value)) +
      geom_boxplot(aes(fill=variable),outlier.shape=NA) +
      theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin,ylimsMax)) +
      annotate(geom=\"segment\",x=0,xend=8,y=ylines,yend=ylines,lty=2) +
      ylab(ylabs) + xlab(\"Samples\") + ggtitle(mytitle2)
   print(p)
   dev.off()
}
";

LOG($outLog, date() . "4. $YW Running R script $LCY$resDir/RESULT.R$N\n");
open (my $outR, ">", "$resDir/RESULT.R") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/RESULT.R$N: $!\n");
print $outR $Rscript;
close $outR;
system("R --vanilla --no-save < $resDir/RESULT.R > $resDir/RESULT.R.LOG 2>&1") == 0 or DIELOG($outLog, date() . " Failed to run $LCY$resDir/RESULT.R$N: $!\n");
my $tailR = `tail $resDir/RESULT.R.LOG`;

LOG($outLog, "\n\n" . date() . " ${LGN}SUCCESS on running $LCY$resDir/RESULT.R$YW.\nLast 5 rows of log message:$N\n$N$tailR


${YW}To Run R script manually, do:$N
R --vanilla --no-save < $resDir/RESULT.R
");
}

####### PARAMETERS
sub get_cluster {
	my ($bedFile, $coreOutFile, $cluster, $bedFileCount, $outLog) = @_;
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($clusterFile) = $bedFile;
	$clusterFile =~ s/.clust.bed/.clust/;
	my $maxClust = -1;
	LOG($outLog, date() . "  1.$bedFileCount. ${LGN}EXIST  $N: clusterFile=$LCY$clusterFile$N\n") if -e $clusterFile;
	LOG($outLog, date() . "  1.$bedFileCount. ${LRD}MISSING$N: clusterFile=$LCY$clusterFile$N\n") if not -e $clusterFile;
	if (-e $clusterFile) {
		my $linecount = -1;
		my $printMOD = 0; #to print example if there are 1 read with 2 clusters.
		LOG($outLog, date() . "   coreOutFile = $YW$coreOutFile$N\n");
		open (my $clusterIn, "<", $clusterFile) or DIELOG($outLog, "Failed to read from clusterFile $clusterFile: $!\n");
		while (my $line = <$clusterIn>) {
			chomp($line);
			$linecount ++;
			my ($currid, $x, $xmax, $y, $ymax, $clust, $curridlong) = split("\t", $line);
			$maxClust = $clust if $maxClust < $clust;
			my $id = $currid;
			LOG($outLog, date() . "NEW id=$id clust=$clust x=$x, xmax=$xmax, y=$y, ymax=$ymax, clust=$clust, $coreOutFile,$id,$clust\n",$DEBUG);# if $id eq "1704132300151533640";
			$cluster->{$coreOutFile}{$id}{$x}{$xmax} = $clust;
		}
	}
#commentcut
	LOG($outLog, date() . "$LCY$bedFile$N cluster=$LGN$maxClust$N\n",$DEBUG);
	#print "temp=$LPR$coreOutFile$N\n";
	return($cluster);
}

sub preprocess_bed {
	my ($bedFile, $coreOutFile, $cluster, $outLog) = @_; 

	my %clust;
	foreach my $id (sort keys %{$cluster->{$coreOutFile}}) {
		foreach my $x (sort {$a <=> $b} keys %{$cluster->{$coreOutFile}{$id}}) {
			foreach my $xmax (sort {$a <=> $b} keys %{$cluster->{$coreOutFile}{$id}{$x}}) {
				my $cluster = $cluster->{$coreOutFile}{$id}{$x}{$xmax};
				my $peakind = defined $clust{$id} ? scalar(@{$clust{$id}}) : 0;
				push(@{$clust{$id}}, $cluster);
				print "peakind = $peakind, id=$id, clust = $cluster, clustid = $clust{$id}[$peakind]\n" if $id eq "160130030742360360" or $id eq "1601300307421077340";
			}
		}
	}
	LOG($outLog, date() . " Getting fasta and calculating GC skew from: $LCY$bedFile$N\n", $DEBUG);
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my $window = 100;
	my $window2 = 200;

	my $bedFileChanged = "$resDir/.TEMP/$bedFilename.BED";

	my %bed;
	open (my $bedIn, "<", $bedFile) or DIELOG($outLog, date() . "Failed to read from $bedFile: $!\n");
	open (my $bedOut, ">", $bedFileChanged) or DIELOG($outLog, date() . "Failed to read from $bedFileChanged: $!\n");
	while (my $line = <$bedIn>) {
		chomp($line);
		my ($chr, $beg, $end, $name, $zero, $strand) = split("\t", $line);
		my ($gene, $id) = $name =~ /^(.+)\.(.+)$/;
		DIELOG($outLog, date() . "Failed to parse gene and id name form bedfile=$bedFile\n") if not defined $gene or not defined $id;
		#my ($id, $id1, $id2, $id3) = parse_readName($id, $outLog);
		#$id = "$id1$id2$id3";
		$bed{$id}{coor}{$beg}{$end} = "$chr\t$beg\t$end\t$gene.$id.MYCLUSTER\t$zero\t$strand";
	}
	close $bedIn;
	foreach my $id (sort keys %bed) {
		my $peakind = -1;
		foreach my $beg (sort {$a <=> $b} keys %{$bed{$id}{coor}}) {
			foreach my $end (sort {$a <=> $b} keys %{$bed{$id}{coor}{$beg}}) {
				$peakind ++;
				my $line = $bed{$id}{coor}{$beg}{$end};
				my $clust = $clust{$id}[$peakind];
				#print "id=$id ind=$peakind clust=$clust\n";# if $id eq "160130030742360360";
				DIELOG($outLog, date() . "Failed to get cluster on bedfile=$bedFile, id=$id beg=$beg end=$end peakind=$peakind\n") if not defined $clust;
				$line =~ s/MYCLUSTER\t/$clust\t/;
				print $bedOut "$line\n";
			}
		}
	}
	close $bedOut;
	my $outputA = "$resDir/.TEMP/$bedFilename\_$window\_A.temp";
	my $outputB = "$resDir/.TEMP/$bedFilename\_$window\_B.temp";
	my $outputC = "$resDir/.TEMP/$bedFilename\_$window\_C.temp";
	my $outputD = "$resDir/.TEMP/$bedFilename\_$window\_D.temp";
	my $outputE = "$resDir/.TEMP/$bedFilename\_$window\_E.temp";
	my $outputF = "$resDir/.TEMP/$bedFilename\_$window\_F.temp";
	my $outputW = "$resDir/.TEMP/$bedFilename\_$window\_W.temp";

	print "$footLoopScriptsFolder/lib/bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFileChanged -o $outputA > $outputA.LOG 2>&1\n";
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFileChanged -o $outputA > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -a -x -$window -y $window -i $bedFileChanged -o $outputB > $outputB.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -a -x 0 -y $window2 -i $bedFileChanged -o $outputC > $outputC.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -x 0 -y 0 -i $bedFileChanged -o $outputW > $outputW.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -b -x -$window2 -y 0 -i $bedFileChanged -o $outputD > $outputD.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -b -x -$window -y $window -i $bedFileChanged -o $outputE > $outputE.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
	system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -b -x 0 -y $window2 -i $bedFileChanged -o $outputF > $outputF.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run $footLoopScriptsFolder/lib/bedtools_bed_change.pl: $!\n");
}




__END__
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
__END__
#			if (not defined $cluster->{$coreOutFile}{$id}) {
#			}
#			elsif (defined $cluster->{$coreOutFile}{$id}) {
#				next if $cluster->{$coreOutFile}{$id}{len} >= $xmax - $x;
#				my $currlen = $xmax - $x;
#				LOG($outLog, date() . "MOD id=$id clust=$cluster->{$coreOutFile}{$id}{clust}, newclust=$clust, num=$number len=$cluster->{$coreOutFile}{$id}{len} < $currlen\n",$DEBUG);
#				LOG($outLog, date() . "Example MOD:\n" . date() . "MOD id=$id clust=$cluster->{$coreOutFile}{$id}{clust}, newclust=$clust, num=$number len=$cluster->{$coreOutFile}{$id}{len} < $currlen\n") if $printMOD == 0;
#				$printMOD = 1;
#				$cluster->{$coreOutFile}{$id}{'$x,$xmax'}{clust} = $clust;
#				$cluster->{$coreOutFile}{$id}{'$x,$xmax'}{len} = ($xmax - $x);
#			}
		}
	}

__END__
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

