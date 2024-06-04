#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_w $opt_g $opt_G $opt_v $opt_n $opt_r $opt_R $opt_B $opt_c $opt_F $opt_0 $opt_J);
getopts("n:vg:w:G:r:R:B:cF0J:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
}

use myFootLib;
use FAlite;

my $homedir = $ENV{"HOME"};

my ($footLoop_script_folder, $version, $md5script) = check_software();
my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
$version_small = $version if not defined $version_small;

my ($max_parallel_run) = defined $opt_J ? $opt_J : 1;

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N -n$LCY <footPeak output directory>$N

${LGN}Optionals$N:

-G $LGN<gene to process>$N] 

-c: Include Cytosine from CpG (use this for in vitro data)


-r Option to run R scripts $LCY(PNG)$N
   -r 0: do not run any R scripts
   -r 1: run only relevant R scripts (default)
   -r 2: run ALL R scripts

-R Option to run R scripts $LCY(PDF)$N
   -R 0: do not run any R scripts (default)
   -R 1: run only relevant R scripts
   -R 2: run ALL R scripts

-B: <BED3 file> add option to add box in the graph

";

die $usage if not defined $opt_n;
die "\nERROR: -n footPeak dir $LCY$opt_n$N doesn't exists!\n\nUsage: $YW$0$N -n <footPeak output directory>\n\n" if not -d $opt_n;

my ($currMainFolder) = `pwd`; chomp($currMainFolder);
$opt_r = 1 if not defined $opt_r;
$opt_R = 0 if not defined $opt_R;

die "\nERROR: -r has to be 0, 1, or 2! (current: $opt_r)\n\n" if defined $opt_r and $opt_r !~ /^[012]$/;

my ($user) = $homedir =~ /home\/(\w+)/;
$user = "USER" if not defined $user;
my $uuid = getuuid();
my $date = date();

my $boxFile = "";
if (defined $opt_B and -e $opt_B) {
	$boxFile = getFullpath($opt_B);
}

main($opt_n);

sub main {
	my @Rscript;
	my ($resDir) = @_;
	my $resDir2 = getFullpath($resDir);
	my ($resDirFullpath) = getFullpath($resDir);
	my %scp;
   makedir("$resDir/.CALL") if not -d "$resDir/\.CALL";
   makedir("$resDir/PEAKS_GENOME") if not -d "$resDir/PEAKS_GENOME";
   makedir("$resDir/PEAKS_LOCAL") if not -d "$resDir/PEAKS_LOCAL";
   makedir("$resDir/PNG") if not -d "$resDir/PNG";
	
	my ($OUTDIRS, $PEAK, $TEMP, $RCONV, $CPG, $ALL);
	($OUTDIRS->{PNG}, $PEAK, $TEMP, $RCONV, $CPG, $ALL) = makeOutDir($resDirFullpath . "/PNG/");
	($OUTDIRS->{PDF}) = makeOutDir($resDirFullpath . "/PDF/");
	my %files;
	my ($footPeak_logFile) = "$resDir/footPeak_logFile.txt";
	open (my $outLog, ">", "$resDir/footPeak_graph_logFile.txt") or die "\n\nFailed to write to $resDir/footPeak_graph_logFile.txt: $!\n\n";

	############
	# LOG START
	LOG($outLog, ">footPeak_graph.pl version $version\n");
	LOG($outLog, ">UUID: $uuid\n", "NA");
	LOG($outLog, ">Date: $date\n", "NA");
	LOG($outLog, ">Run script: $0 -n $opt_n\n", "NA");
	
	############

	my %coor;
	my @lines   = `cat $footPeak_logFile`;
	LOG($outLog, " Parsing footpeak logfile $footPeak_logFile\n");
	my ($label) = `cat $resDir/.LABEL`; chomp($label);
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_logFile!\n\n") if not -e $footPeak_logFile;
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $resDir/.LABEL!\n\n") if not -e "$resDir/.LABEL";
	my ($thres, $window);

LOG($outLog, "

If R dies for any reason, make sure you have these required R libraries:
- RColorBrewer v1.1-2
- gridExtra v2.3
- labeling v0.3
- reshape2 v1.4.3
- ggplot2 v3.1.0
- GMD v0.3.3
");

	foreach my $line (@lines) {
		chomp($line);
		if ($line =~ /^[ \t]*def=.+, coor=.+/) {
			$line =~ s/(^\s+|\s+$)//g;
			my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
			LOG($outLog, "gene=$gene,chr=$CHR,beg=$BEG,end=$END,gene=$GENE,val=$VAL,strand=$STRAND\n","NA");
	   	if (defined $opt_G and $gene !~ /$opt_G/i) {
	   	   LOG($outLog, date() . " Skipped $LCY$gene$N as it doesn't contain $LGN-G $opt_G$N\n");
	   	   next;
	   	}
			$GENE = uc($GENE);
			DIELOG($outLog, "\n\ndied at processing $LCY$footPeak_logFile$N: can't parse index file def gene lqines\n\n$line\n\n") if not defined $STRAND;
			%{$coor{$GENE}} = ("CHR" => $CHR, "BEG" => $BEG, "END" => $END, "VAL" => $VAL);
			$coor{$GENE}{STRAND} = $STRAND eq "+" ? "Pos" : $STRAND eq "-" ? "Neg" : $STRAND =~ /^(Pos|Neg|Unk)$/ ? $STRAND : "Unk";
		}
		elsif ($line =~ /^-t thrshld\s+:/) {
			($thres) = $line =~ /^-t thrshld\s+:\s+(\-?\d+\.?\d*)$/;
			$thres = "0." . $thres if $thres > 1;
		}
		elsif ($line =~ /^-w window\s+:/) {
			($window) = $line =~ /^-w window\s+:\s+(\-?\d+\.?\d*)$/;
		}
	}
	my %gene;
	my @files;
	my @types = qw(CH CG GH GC);
	foreach my $GENE (sort keys %coor) {
		my $mygene = $GENE;
		my @strands = qw(Pos Neg);
		for (my $h1 = 0; $h1 < @strands; $h1++) {
			for (my $h2 = 0; $h2 < 4; $h2++) {
				my $strand = $strands[$h1];
				my $type = $types[$h2];
				my $thres2 = $thres; $thres =~ s/0$//;
				my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK";
				my $nopkFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.NOPK";
				#LOG($outLog, "$resDir/PNG/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK.out.PEAK.png.ALL.conly.png\n");
				#LOG($outLog, "$resDir/PNG/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.NOPK.out.NOPK.png.ALL.conly.png\n");
				my $peakflag = getFlag($peakFile, "Pos", $strand, $type);
				my $nopkflag = getFlag($nopkFile, "Pos", $strand, $type);
				if ($peakflag !~ /RCONV/) {
					LOG($outLog, "${LPR}PNGOUT $peakflag:$N\n");
					LOG($outLog, "scp mitochi\@franklin.hpc.ucdavis.edu:/$resDir2/PNG/$peakflag/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK.out.$peakflag.png.ALL.conly.png ./\n");
					LOG($outLog, "scp mitochi\@franklin.hpc.ucdavis.edu:/$resDir2/PNG/$nopkflag/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.NOPK.out.$nopkflag.png.ALL.conly.png ./\n");
				}
				#next if $peakFile !~ /PBEH2_BCBC4_PLASMIDPFC9NTBSPQI13_DESCTXLINEARBSAINTBSPQI13/;
				$peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres2\_$type.PEAK" if not -e $peakFile;
				die "Can't find peakFile $LCY$peakFile$N\n" if not -e $peakFile;
				LOG($outLog, "label=$label, gene$mygene, strand=$strand, peak=$peakFile\n","NA");
				$files{$peakFile} = $mygene;
				if ($opt_r eq 0 or $opt_r eq 1) {
					if (defined $opt_c) {
						next if $peakFile !~ /(Pos.+CG|Pos.+CH|Neg.+GC|Neg.+GH)/;
					}
					else {
						next if $peakFile !~ /(Pos.+CH|Neg.+GH)/;
					}
					push(@files, $peakFile);
				}
				elsif ($opt_r eq 2) {
					push(@files, $peakFile);
				}
				#last if @files > 5;
			}
			#last if @files > 5;
		}
		#last if @files > 5;
	}
	my %Rscripts; 
	my $lastfile = -1; #debug
	DIELOG($outLog, "\n\nERROR: There is no files in $LCY\%files$N defined!\n") if (keys %files) == 0;
	my $fileCount = 0;
	my $totalFile = (keys %files);
	my $lastGENE = -1;

	my $sbatch_these_cmd = "footPeak_graph_sbatch_2.pl -n $resDir -i FILENAME -r $opt_r -R $opt_R";
	$sbatch_these_cmd .= " -c" if defined $opt_c;
	$sbatch_these_cmd .= " -0" if defined $opt_0;
	$sbatch_these_cmd .= " -F" if defined $opt_F;
	$sbatch_these_cmd .= " -B $boxFile" if defined $opt_B;
	#my $cmd = "$footLoop_script_folder/footPeak_sbatch_2.pl $indexFile $seqFile FNINDICE FILENAME $outDir >
	my $force_sbatch = 1 if defined $opt_F;
	my $outsbatchDir = "$resDir/.footPeak_graph_sbatch/";
	system("mkdir -p $outsbatchDir") if not -d $outsbatchDir;
	
	sbatch_these($sbatch_these_cmd, "footPeak_graph_sbatch", \@files, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir);
}


###############
# Subroutines #
###############

sub sbatch_these {
   #my ($cmd, $suffix, $ext, $filesARRAY, $max_parallel_run, $outLog, $force_sbatch, $folderwant) = @_;
   my ($cmd, $suffix, $filesARRAY, $max_parallel_run, $outLog, $force_sbatch, $folderwant) = @_;

   my %force;

   # - suffix : $file\_$suffix.sbatch/sbout/folder
   # - ext    : determine donefile (<$file\_$suffix/*.$ext>)

   my $jobidhash;
   my $totalfiles = scalar(@{$filesARRAY});
   for (my $i = 0; $i < @{$filesARRAY}; $i++) {

      my $file = $filesARRAY->[$i];
      my ($folder, $filename) = getFilename($file, "folderfull");
      $folderwant = $folder if not defined $folderwant;
      my $sbatchfile = "$folderwant/$filename\_$suffix.sbatch";
      my $sboutfile  = "$folderwant/$filename\_$suffix.sbout";
      my $donefile   = "$folderwant/$filename\_$suffix.done";

      my $cmdcopy = $cmd;
         $cmdcopy =~ s/FILENAME/$file/g;

      if ($cmdcopy !~ /FOLDER/) {
         system("mkdir -p $folderwant/$filename\_$suffix") if not -e "$folderwant/$filename\_$suffix";
      }
      else {
         $cmdcopy =~ s/FOLDER/$folderwant/g;
      }
      if ($i == 0) {
         print_cmd($cmdcopy, $outLog);
      }

      if ($cmd =~ /FNINDICE/) {
         $cmdcopy =~ s/FNINDICE/$i/g;
      }

      $sboutfile =~ s/\/+/\//g;
      my $sbatchprint = "";
         $sbatchprint .= "#!/bin/bash -l\n";
         $sbatchprint .= "#SBATCH -n 2 -N 1 -p high --mem 4000 -t 999:99:99\n";
         $sbatchprint .= "#SBATCH --job-name \"$filename\_$suffix\"\n";
         $sbatchprint .= "#SBATCH --output \"$sboutfile\"\n\n";
         $sbatchprint .= "conda activate footLoop2\n";
         $sbatchprint .= "$cmdcopy && echo \"Done!\" > $donefile\n\n";

      #"sbatchFile=\n$LCY$sbatchfile$N\n\n" if defined $opt_0;
      open (my $out, ">", $sbatchfile) or die "Can't write to $LCY$sbatchfile$N: $!\n";
      print $out $sbatchprint;
      close $out;

      my $iprint = $i + 1;
      if (not defined $force{0} and not defined $force_sbatch and -e $donefile) {
         LOG($outLog, "\n" . date() . "${LPR}$iprint/$totalfiles sbatch_these $suffix$N: sbatch $LCY$sbatchfile$N # ${LGN}DONE$N\n");
         next;
      }
      else {
         LOG($outLog, "\n" . date() . "${LPR}$iprint/$totalfiles sbatch_these $suffix$N: sbatch: $LCY$sbatchfile$N\n");
      }

      if (defined $opt_0) { # Debug
         next;
      }

      if ($i != 0) {
         my $sleep = 0;
         while (1) {
            last if $i < $max_parallel_run;
            my ($job_left) = squeue_check($jobidhash);
            LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
            last if ($job_left < $max_parallel_run);
            $sleep ++;
            sleep 5;
         }
      }
      my ($jobid) = `sbatch $sbatchfile`;
      chomp($jobid);
      ($jobid) = $jobid =~ /^Submi.+job (\d+)$/;
      next if not defined $jobid;
      $jobidhash->{$jobid} = 1;
      LOG($outLog, "$YW$i$N $LCY$filename$N $LGN$jobid$N\n");
      #system("touch $file.done") == 0 or die "failed to touch $file.done: $!\n";
   }
   my $sleep = 0;
   #while (1) {
   #  my ($job_left) = squeue_check($jobidhash);
   #  LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 60 == 0;
   #  last if ($job_left < $max_parallel_run);
   #  $sleep ++;
   #  sleep 1;
   #}
   #$sleep = 0;
   while (1) {
      my ($job_left) = squeue_check($jobidhash);
      LOG($outLog, "\n" . date() . "$job_left jobs left!\n") if $sleep % 12 == 0;
      last if $job_left == 0;
      $sleep ++;
      sleep 5;
   }
   LOG($outLog, "\n" . date() . "All have been run!\n\n");
   return(0);
}

sub squeue_check {
   my ($jobidhash, $outLog) = @_;
   my @squeue = `squeue`;
   my $squeuehash;
   foreach my $line (@squeue) {
      next if $line =~ /JOBID\s+PARTITION.+/;
      my ($jobid) = $line =~ /^\s*(\d+)\s+/;
      if (not defined $jobid) {
         LOG($outLog, "Can't parse jobid from line=$LCY$line$N\n");
         next; # just next so we don't kill the script...
      }
      next if not defined $jobidhash->{$jobid};
      $squeuehash->{$jobid} = 1;
   }
   foreach my $jobid (keys %{$jobidhash}) {
      next if defined $squeuehash->{$jobid};
      undef $jobidhash->{$jobid};
      delete $jobidhash->{$jobid};
   }
   my ($total) = scalar(keys %{$jobidhash});
   return ($total);
}

sub print_cmd {
   my ($cmd, $outLog) = @_;
   #LOG($outBigCMD, "\n$cmd\n","NA");
   if ($cmd !~ /^#/) {
      $cmd =~ s/^/    /;
      $cmd =~ s/ \-/ \\\n      \-/g;
   }
   else {
      $cmd =~ s/^#/   # /;
      $cmd =~ s/ \-/ \\\n     # \-/g;
   }
   LOG($outLog, "$LGN\n$cmd\n$N\n");
   #LOG($outBigCMD, "\n$cmd\n","NA");
}


sub Rscript {
	my ($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout, $resDir) = @_;
	my $currFileID = $currFile . ".id";
	my $R;
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
	($bedFilename) =~ s/(CG|CH|GC|GH).+$/$1/;
	my $peakminVal = $currFile =~ /PEAK.out$/ ? 8 : 6;
	my $plotminReads = 4;

# -------------------- $R->{readTable}
	$R->{readTable} = "	

#plotminReads = $plotminReads

theme_blank <- theme_bw()
theme_blank\$line <- element_blank()
theme_blank\$rect <- element_blank()
theme_blank\$strip.text <- element_blank()
theme_blank\$axis.text <- element_blank()
theme_blank\$axis.ticks <- element_blank()
theme_blank\$axis.title <- element_blank()
theme_blank\$legend.title=element_blank()
theme_blank\$axis.line = element_line()
theme_blank\$plot.margin = structure(c(0, 0, 0,0), unit = \"lines\", valid.unit = 3L, class = \"unit\")
theme_blank\$panel.grid=element_blank()
theme_blank\$panel.grid.major=element_blank()
theme_blank\$panel.grid.minor=element_blank()
theme_blank\$panel.background=element_blank()
theme_blank\$panel.border= element_blank()
theme_blank\$legend.position=\"none\"

p.png.scale = 1
p.pdf.scale = 0.33

#####################
# Read Table
df = read.table(\"$currFile\",sep=\"\\t\")
df.id = read.table(\"$currFileID\",sep=\"\\t\",colClasses=c(\"factor\",\"factor\"))
print(head(df[1:10]))
print(head(df.id))
colnames(df) = c(\"V1\",seq(1,dim(df)[2]-1))
colnames(df.id) = c(\"V1\",\"ID\")
df\$id = df.id\$ID
cluster_color = c()

#####################
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

";

# -------------------- $R->{clusterFile}

	$R->{clusterFile} = "

#####################
# Cluster
clust = read.table(\"$curr_cluster_file\",header=F,sep=\"\\t\",colClasses=c(\"factor\",\"integer\",\"integer\",\"integer\",\"integer\",\"numeric\",\"character\"))
colnames(clust) = c(\"id\",\"x\",\"xmax\",\"y\",\"ymax\",\"clust\",\"id2\")
print(head(clust))
clust\$y = seq(1,dim(clust)[1])
clust2 = as.data.frame(aggregate(clust\$y,by=list(clust\$clust,clust\$id2),min))
clust2\$ymax = aggregate(clust\$y,by=list(clust\$clust,clust\$id2),max)\$x
clust2\$xpos0 = aggregate(clust\$x,by=list(clust\$clust,clust\$id2),min)\$x
clust2\$xpos1 = aggregate(clust\$xmax,by=list(clust\$clust,clust\$id2),max)\$x
colnames(clust2) = c(\"clust\",\"id2\",\"ymin\",\"ymax\",\"xpos0\",\"xpos1\")
clust2 = subset(clust2,select=c(\"clust\",\"ymin\",\"ymax\",\"xpos0\",\"xpos1\",\"id2\"))
print(clust2)
clust2\$xmin = 1
clust2\$xmax = 140
clust2\$clust = clust2\$clust + 10
clust = subset(clust,select=-id2)
clust = subset(clust,select=c(\"id\",\"y\",\"clust\"))
print(head(clust))
print(head(clust2))
df3 = merge(df,clust,by=\"id\")
print(head(df3[1:10]))
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
		} else if (peak != 0 & peak/2 >= conv) {
			df4[x,y] = as.integer(length(a[a == 8 | a == 9]) / length(a) * 9+0.5)
		} else if (conv != 0 & conv > peak/2) {
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
colnames(df3) = c(\"y\",\"clust\",\"x\",\"value\")
df3\$x     = as.numeric(as.character(df3\$x))
df3\$y     = as.numeric(as.character(df3\$y))
df3\$clust = as.numeric(as.character(df3\$clust))
df3\$value = as.numeric(as.character(df3\$value))
df5 = df3; df5[df5\$value > 10 | df5\$value < -10,]\$value = 0 #turn non-peak heatmap value into 0
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
clust = subset(clust,select=c(\"id\",\"y\"))
df = merge(df,clust,by=\"id\")
df = subset(df,select=-id)
df3\$x = as.numeric(as.character(df3\$x))
df3\$y = as.numeric(as.character(df3\$y))
df3\$value  = as.numeric(as.character(df3\$value))
df4\$clust2 = df4\$clust + 10
print(unique(df4\$clust2))
greens = rev(brewer.pal(9,\"Greens\"))
reds = brewer.pal(9,\"Reds\")
p3.col=c(
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

get_cluster_color = function(cluster) {
	#9 SET1 color
	cluster_color_array = c(\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\",\"#ffff33\",\"#a65628\",\"#f781bf\",\"#999999\")
	cluster_color_array = c(cluster_color_array,cluster_color_array,cluster_color_array)
   cluster_color_values = c()
   cluster_color_names = c()
   for (i in 1:length(unique(cluster))) {
      curr.clust = cluster[i];
      cluster_color_array.i = i \%\% 9;
      if (cluster_color_array.i == 0) {
         cluster_color_array.i = 9
      }
      cluster_color_names[i] = cluster_color_array[cluster_color_array.i]
      cluster_color_values[i] = curr.clust
   }
   cluster_color = cluster_color_names
   names(cluster_color)= cluster_color_values
   cluster_color
}
cluster_color = get_cluster_color(unique(df4\$clust2))
cluster_length = length(unique(df4\$clust2))
print(cluster_color)

p3.png.scale = p.png.scale
p3.pdf.scale = p.pdf.scale

p3.png = ggplot(df3,aes(x,y)) +
	geom_tile(aes(fill=as.factor(value))) +
	geom_rect(data=df5,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x=xmin,y=ymin),
				 size=0.8*p3.png.scale,fill=rgb(0,0,0,alpha=0),color=rgb(0,0,0,alpha=0.25) ) + # SIZE
	geom_rect(data=df4,aes(xmin=xmin,xmax=xmax,ymin=ymin,
		ymax=ymax,x=xmin,y=ymin,color=as.factor(df4\$clust2)),
		size=1*p3.png.scale,fill=rgb(0,0,0,alpha=0)) + # SIZE
	geom_rect(data=df4,aes(fill=as.factor(df4\$clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax),size=0.5*p3.png.scale) + # SIZE
	geom_text(data=df4,aes(x=10,y=(ymin+ymax)/2,label=clust2-10),hjust=0,size=5*p3.png.scale) + # SIZE
	theme_bw() + theme(legend.position=\"none\") + coord_cartesian(ylim=c(-1,cluster_length+2)) +
	scale_fill_manual(values=c(p3.col,cluster_color)) +
	scale_color_manual(values=c(p3.col,cluster_color)) +
	scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	theme_blank


p3.pdf = ggplot(df3,aes(x,y)) +
	geom_tile(aes(fill=as.factor(value))) +
	geom_rect(data=df5,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x=xmin,y=ymin),
				 size=0.8*p3.pdf.scale,fill=rgb(0,0,0,alpha=0),color=rgb(0,0,0,alpha=0.25) ) + # SIZE
	geom_rect(data=df4,aes(xmin=xmin,xmax=xmax,ymin=ymin,
		ymax=ymax,x=xmin,y=ymin,color=as.factor(df4\$clust2)),
		size=1*p3.pdf.scale,fill=rgb(0,0,0,alpha=0)) + # SIZE
	geom_rect(data=df4,aes(fill=as.factor(df4\$clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax),size=0.5*p3.pdf.scale) + # SIZE
	geom_text(data=df4,aes(x=10,y=(ymin+ymax)/2,label=clust2-10),hjust=0,size=3) + # SIZE
	theme_bw() + theme(legend.position=\"none\") + coord_cartesian(ylim=c(-1,cluster_length+2)) +
	scale_fill_manual(values=c(p3.col,cluster_color)) +
	scale_color_manual(values=c(p3.col,cluster_color)) +
	scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	theme_blank

";

# -------------------- $R->{noclusterFile}
	$R->{noclusterFile} = " 

##########################
# Cluster: NO CLUSTER
df\$y = seq(1,dim(df)[1])\ndf=subset(df,select=-id)

";

# -------------------- $R->{peakbedFile}

	$R->{peakbedFile} = "

##########################
# Peak BedFile
bed = read.table(\"$bedFile\",sep=\"\\t\")
bed = merge(subset(df,select=c(\"V1\",\"y\")),bed,by=\"V1\")

";

# -------------------- $R->{mainplot}
	$R->{mainplot} .= "
##########################
# Main Plot Part 1
dm = melt(df,id.vars=c(\"V1\",\"y\"))
dm\$variable = as.numeric(as.character(dm\$variable))
p.col = c(
			\"0\"=\"#f0f0f0\",
			\"1\"=\"white\",
			\"4\"=\"cornsilk\",
			\"5\"=\"cornsilk\",
			\"6\"=\"green4\",
			\"7\"=\"seagreen4\",
			\"8\"=\"red4\",
			\"9\"=\"maroon4\"
)
if (length(cluster_color) > 0) {
	p.col = c(p.col, cluster_color)
}
mywindow = 1000
p = ggplot(dm,aes(variable,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL) + theme_blank +
	 ggtitle(paste(\"(peak=\",$totpeak,\"; nopk=\",$totnopk,\")\",sep=\"\"))
p.heatmaponly = ggplot(dm,aes(variable,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL) + theme_blank
p.png = p
p.pdf = p
";

	$R->{box} = "

	box = read.table(\"$boxFile\",sep=\"\\t\")
	p.png = p.png + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.pdf = p.pdf + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	
	";

	$R->{box_nopk} = "

	box = read.table(\"$boxFile\",sep=\"\\t\")
	p.rand.100.png  = p.rand.100.png  + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.100\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.rand.100.pdf  = p.rand.100.pdf  + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.100\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	p.rand.1000.png = p.rand.1000.png + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.rand.1000.pdf = p.rand.1000.pdf + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	
	";

# -------------------- $R->{mainplot_nopk_rand_100}
	$R->{mainplot_nopk_rand_100} .= "

##########################
# Main Plot Part 1 nopk subset 100

# Randomly take reads
df.total = dim(df)[1]
set.seed(42)
if (df.total > 100) {
	df.rand.100 = df[sample( seq(1,df.total) ,100,replace=F),]
	print(\"randoming 100. DF TOTAL = \"); print(df.total)
} else {
	df.rand.100 = df
	print(\"NOT randoming 100. DF TOTAL = \"); print(df.total)
}

# sorting by hclust
if (dim(df.rand.100)[1] < 1000) {
	h = hclust(dist(df.rand.100[,-1]))
	df.rand.100 = df.rand.100[h\$order,]
} else if (dim(df.rand.100)[2] < 10) {
	mysum = apply(df.rand.100[,-1],1,sum)
	df.rand.100 = df.rand.100[order(mysum),]
}
df.rand.100\$y = seq(1,dim(df.rand.100)[1])

write.table(df.rand.100,file=\"$currFile.100.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
print(\"Wrote to $currFile.100.rand\")

dm.rand.100 = melt(df.rand.100,id.vars=c(\"V1\",\"y\"))
print(dim(df.rand.100))
print(dim(dm.rand.100))
dm.rand.100\$variable = as.numeric(as.character(dm.rand.100\$variable))

p.rand.100 = ggplot(dm.rand.100,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 theme(
	 	line = element_blank(),
	 	axis.text = element_blank(),
	 	axis.title = element_blank()
	 ) + 
	 ggtitle(paste(\"(peak=\",$totpeak,\"; nopk=\",$totnopk,\")\",sep=\"\"))

p.heatmaponly.rand.100 = ggplot(dm.rand.100,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 theme(
	 	line = element_blank(),
	 	axis.text = element_blank(),
	 	axis.title = element_blank()
	 )

	p.rand.100.pdf = p.rand.100 + theme(plot.title = element_text(size = 10*p.pdf.scale))
	p.rand.100.png = p.rand.100 + theme(plot.title = element_text(size = 10*p.png.scale))
";


# -------------------- $R->{mainplot_nopk_rand_1000}
	$R->{mainplot_nopk_rand_1000} .= "

##########################
# Main Plot Part 1 nopk rand
df.total = dim(df)[1]
set.seed(42)
if (df.total > 1000) {
	df.rand.1000 = df[sample( seq(1,df.total) ,1000,replace=F),]
	print(\"randoming 1000. DF TOTAL = \"); print(df.total)
} else {
	df.rand.1000 = df
	print(\"NOT randoming 1000. DF TOTAL = \"); print(df.total)
}

# sorting by hclust
if (dim(df.rand.1000)[1] < 1000) {
	h = hclust(dist(df.rand.1000[,-1]))
	df.rand.1000 = df.rand.1000[h\$order,]
} else if (dim(df.rand.1000)[2] < 10) {
	mysum = apply(df.rand.1000[,-1],1,sum)
	df.rand.1000 = df.rand.1000[order(mysum),]
}
df.rand.1000\$y = seq(1,dim(df.rand.1000)[1])

write.table(df.rand.1000,file=\"$currFile.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
print(\"Wrote to $currFile.rand\")

dm.rand.1000 = melt(df.rand.1000,id.vars=c(\"V1\",\"y\"))
print(dim(df.rand.1000))
print(dim(dm.rand.1000))
dm.rand.1000\$variable = as.numeric(as.character(dm.rand.1000\$variable))

p.rand.1000 = ggplot(dm.rand.1000,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 theme(
	 	line = element_blank(),
	 	axis.text = element_blank(),
	 	axis.title = element_blank()
	 ) + 
	 ggtitle(paste(\"(peak=\",$totpeak,\"; nopk=\",$totnopk,\")\",sep=\"\"))

p.heatmaponly.rand.1000 = ggplot(dm.rand.1000,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 theme(
	 	line = element_blank(),
	 	axis.text = element_blank(),
	 	axis.title = element_blank()
	 )

	p.rand.1000.pdf = p.rand.1000 + theme(plot.title = element_text(size = 10*p.pdf.scale))
	p.rand.1000.png = p.rand.1000 + theme(plot.title = element_text(size = 10*p.png.scale))

";

# -------------------- $R->{mainplot_nopk}
	$R->{mainplot_nopk} .= "

##########################
# Main Plot Part 1 nopk

dm = melt(df,id.vars=c(\"V1\",\"y\"))
dm\$variable = as.numeric(as.character(dm\$variable))
p.col = c(
			\"0\"=\"#f0f0f0\",
			\"1\"=\"white\",
			\"4\"=\"cornsilk\",
			\"5\"=\"cornsilk\",
			\"6\"=\"green4\",
			\"7\"=\"seagreen4\",
			\"8\"=\"red4\",
			\"9\"=\"maroon4\"
)
if (length(cluster_color) > 0) {
	p.col = c(p.col, cluster_color)
}

mywindow = 1000
plot_list = list();


for (i in seq(1,as.integer(dim(df)[1] / mywindow) + 1)) {
   beg = (i-1) * mywindow + 1
   end = as.integer(dim(df)[1] / mywindow) + 1
   if (i == end) {
      end = beg + dim(df)[1] \%\% mywindow
   } else {
      end = i * mywindow
   }
   dm.temp = dm[dm\$y >= beg & dm\$y <= end,]
   dm.temp\$y = dm.temp\$y - beg + 1
   print(head(dm.temp))
	p = ggplot(dm.temp,aes(variable,y)) +  
		 geom_tile(aes(fill=as.factor(value))) + 
		 theme_bw() + theme(legend.position=\"none\") + 
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0)) +
		 theme(
		 	line = element_blank(),
		 	axis.text = element_blank(),
		 	axis.title = element_blank()
		 ) + 
		 ggtitle(paste(\"PLOT #\",i,\" (total peak=\",$totpeak,\"; total nopk=\",$totnopk,\")\",sep=\"\"))
	plot_list[[i]] = p
	p.png = p
	p.pdf = p
}
";


# -------------------- $R->{mainplotClusterAddition}
	$R->{mainplotClusterAddition} = "

# Main Plot Cluster Addition
# geom rect default size = 0.5
# geom text default size = 10
p.png = p.png + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.png.scale) + # SIZE
	 geom_rect(data=clust2,aes(color=as.factor(clust),x=xpos0,y=ymin,xmin=xpos0,xmax=xpos1,ymin=ymin,ymax=ymax),size=0.5*p.png.scale,fill=rgb(1,1,1,alpha=0),lwd=1*p.png.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=paste(clust-10,\"(\",id2,\")\",sep=\"\")),hjust=0,size=5*p.png.scale) +
	 theme(plot.title = element_text(size = 10*p.png.scale))

p.pdf = p.pdf + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.pdf.scale) + # SIZE
	 geom_rect(data=clust2,aes(color=as.factor(clust),x=xpos0,y=ymin,xmin=xpos0,xmax=xpos1,ymin=ymin,ymax=ymax),size=0.5*p.pdf.scale,fill=rgb(1,1,1,alpha=0),lwd=1*p.pdf.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=paste(clust-10,\"(\",id2,\")\",sep=\"\")),hjust=0,size=10*p.png.scale) +
	 theme(plot.title = element_text(size = 20*p.pdf.scale))

";

# -------------------- 	$R->{mainplotPeakBedAddition}
	$R->{mainplotPeakBedAddition} = "

# Main Plot Peak Bed Addition
if (length(bed) != 0 & dim(bed)[1] > 0) {

	bed\$variable=1
	p.png =	p.png + 
			 	geom_rect(
					data=bed, 
					aes(xmin=V2, xmax=V3, ymin=y-0.5, ymax=y+0.5),
					fill=rgb(0,0,0,alpha=0),
					color=rgb(1,0,0,alpha=0.25),
			 		size=0.5*p.png.scale # SIZE
				)
	p.heatmaponly =	p.heatmaponly + 
			 	geom_rect(
					data=bed, 
					aes(xmin=V2, xmax=V3, ymin=y-0.5, ymax=y+0.5),
					fill=rgb(0,0,0,alpha=0),
					color=rgb(1,0,0,alpha=0.25),
			 		size=0.5*p.png.scale # SIZE
				)
	p.pdf =	p.pdf + 
			 	geom_rect(
					data=bed, 
					aes(xmin=V2, xmax=V3, ymin=y-0.5, ymax=y+0.5),
					fill=rgb(0,0,0,alpha=0),
					color=rgb(1,0,0,alpha=0.25),
			 		size=0.5*p.pdf.scale # SIZE
				)

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
	$R->{secondplotConversionGraph} = "

# Calculate % Conversion
df2 = subset(df,select=c(-V1,-y));
print(head(df))
print(tail(df))
if (length(df2[df2 < $peakminVal]) > 0) {
	df2[df2 < $peakminVal] = 0;
}
if (length(df2[df2 >= $peakminVal]) > 0) {
	df2[df2 >= $peakminVal] = 1
}

min.df2.x2 = min(as.numeric(as.character(df2\$x)))
max.df2.x2 = max(as.numeric(as.character(df2\$x)))
df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean))
df2temp = data.frame(x=NA,x2=NA,y=NA,y2=NA)
#if (dim(df2[df2\$y > 0,])[1] > $plotminReads) {
if (dim(df2[df2\$y > 0,])[1] > $plotminReads & $plotminReads > 5) {
	df2 = df2[df2\$y > 0,]
	df2\$x = as.numeric(as.character(df2\$x));
	df2\$y = as.numeric(as.character(df2\$y));
	df2\$x2 = df2\$x
	df2\$y2 = df2\$y
	for (i in 1:(dim(df2)[1]-10)) {
		a = df2[df2\$x >= df2[i,]\$x & df2\$x <= (df2[i,]\$x+10),]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2[i,]\$x2 = df2[i,]\$x
			df2[i,]\$y2 = mean(a\$y)
			print(paste(\"#\",i))
			print(a\$x)
			print(a\$y)
			print(df2[i,]\$x2)
			print(df2[i,]\$y2)
		}
	}
	mins = seq(1,as.integer(df2[1,]\$x2)-1,10)
	maxs = seq(max(df2\$x2),dim(df)[2]-2,10)
	df2 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2)
	df2 = rbind(df2,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2 = data.frame(x=seq(1,dim(df)[2]), y=0, x2=seq(1,dim(df)[2]), y2=0);
}

df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean))
df2temp = data.frame(x=NA,x2=NA,y=NA,y2=NA)
if (dim(df2[df2\$y > 0,])[1] > $plotminReads & $plotminReads > 5) {
#if (dim(df2[df2\$y > 0,])[1] > $plotminReads) {
	df2 = df2[df2\$y > 0,]
	df2\$x = as.numeric(as.character(df2\$x));
	df2\$y = as.numeric(as.character(df2\$y));
	df2\$x2 = df2\$x
	df2\$y2 = df2\$y
	for (i in 1:(dim(df2)[1]-10)) {
		a = df2[df2\$x >= df2[i,]\$x & df2\$x <= (df2[i,]\$x+10),]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2[i,]\$x2 = df2[i,]\$x
			df2[i,]\$y2 = mean(a\$y)
			print(paste(\"#\",i))
			print(a\$x)
			print(a\$y)
			print(df2[i,]\$x2)
			print(df2[i,]\$y2)
		}
	}
	mins = seq(1,as.integer(df2[1,]\$x2)-1,10)
	maxs = seq(max(df2\$x2),dim(df)[2]-2,10)
	df2 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2)
	df2 = rbind(df2,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2 = data.frame(x=seq(1,dim(df)[2]), y=0, x2=seq(1,dim(df)[2]), y2=0);
}


min.df2.x2 = min(c(df2\$x2))
max.df2.x2 = max(c(df2\$x2))
df2temp = data.frame(x=NA,x2=seq(min.df2.x2,max.df2.x2),y=NA,y2=NA)
for (i in dim(df2temp)[1]) {
	if (dim(df2[df2\$x2 == df2temp[i,]\$x2,])[1] > 0) {
		df2temp[i,]\$x = df2[i,]\$x
		df2temp[i,]\$x2 = df2[i,]\$x2
		df2temp[i,]\$y = df2[i,]\$y
		df2temp[i,]\$y2 = df2[i,]\$y2
	}
}
write.table(df2,file=\"$currFile.fixedwig\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
write.table(df2temp,file=\"$currFile.fixedwig2\",quote=F,row.names=F,col.names=T,sep=\"\\t\")

# P2 % Conversion XY Plot
	p2.png.scale = p.png.scale
	p2.pdf.scale = p.pdf.scale

p2.png = 
	ggplot(df2,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme_blank +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.png.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

p2.pdf = 
	ggplot(df2,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.pdf.scale) +  # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.pdf.scale) + # SIZE 
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme_blank +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.pdf.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{secondplotConversionGraph_rand_100} = "

# Calculate % Conversion
df2.rand.100 = subset(df.rand.100,select=c(-V1,-y));
if (length(df2.rand.100[df2.rand.100 < $peakminVal]) > 0) {
	df2.rand.100[df2.rand.100 < $peakminVal] = 0;
}
if (length(df2.rand.100[df2.rand.100 >= $peakminVal]) > 0) {
	df2.rand.100[df2.rand.100 >= $peakminVal] = 1
}
df2.rand.100 = data.frame(x=seq(1,dim(df2.rand.100)[2]), y=apply(df2.rand.100,2,mean))
if (dim(df2.rand.100[df2.rand.100\$y > 0,])[1] > $plotminReads & $plotminReads > 5) {
#if (dim(df2.rand.100[df2.rand.100\$y > 0,])[1] > $plotminReads) {
	df2.rand.100 = df2.rand.100[df2.rand.100\$y > 0,]
	df2.rand.100\$x = as.numeric(as.character(df2.rand.100\$x));
	df2.rand.100\$y = as.numeric(as.character(df2.rand.100\$y));
	df2.rand.100\$x2 = df2.rand.100\$x
	df2.rand.100\$y2 = df2.rand.100\$y
	for (i in 1:(dim(df2.rand.100)[1]-10)) {
		a = df2.rand.100[df2.rand.100\$x >= df2.rand.100[i,]\$x & df2.rand.100\$x <= df2.rand.100[i+10-1,]\$x,]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2.rand.100[i,]\$y2 = mean(a\$y)
			df2.rand.100[i,]\$x2 = mean(a\$x)
		}
	}
	mins = seq(1,as.integer(df2.rand.100[1,]\$x2)-1,10)
	maxs = seq(max(df2.rand.100\$x2),dim(df.rand.100)[2]-2,10)
	df2.rand.100 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2.rand.100)
	df2.rand.100 = rbind(df2.rand.100,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2.rand.100 = data.frame(x=seq(1,dim(df.rand.100)[2]), y=0, x2=seq(1,dim(df.rand.100)[2]), y2=0);
}
 
# P2 % Conversion XY Plot
p2.rand.png.scale = p.png.scale
p2.rand.pdf.scale = p.pdf.scale

p2.rand.100.png = 
	ggplot(df2.rand.100,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.png.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

p2.rand.100.pdf = 
	ggplot(df2.rand.100,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.pdf.scale) +  # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.pdf.scale) + # SIZE 
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.pdf.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{secondplotConversionGraph_rand_1000} = "

# Calculate % Conversion
df2.rand.1000 = subset(df.rand.1000,select=c(-V1,-y));
if (length(df2.rand.1000[df2.rand.1000 < $peakminVal]) > 0) {
	df2.rand.1000[df2.rand.1000 < $peakminVal] = 0;
}
if (length(df2.rand.1000[df2.rand.1000 >= $peakminVal]) > 0) {
	df2.rand.1000[df2.rand.1000 >= $peakminVal] = 1
}
df2.rand.1000 = data.frame(x=seq(1,dim(df2.rand.1000)[2]), y=apply(df2.rand.1000,2,mean))
if (dim(df2.rand.1000[df2.rand.1000\$y > 0,])[1] > $plotminReads & $plotminReads > 5) {
#if (dim(df2.rand.1000[df2.rand.1000\$y > 0,])[1] > $plotminReads) {
	df2.rand.1000 = df2.rand.1000[df2.rand.1000\$y > 0,]
	df2.rand.1000\$x = as.numeric(as.character(df2.rand.1000\$x));
	df2.rand.1000\$y = as.numeric(as.character(df2.rand.1000\$y));
	df2.rand.1000\$x2 = df2.rand.1000\$x
	df2.rand.1000\$y2 = df2.rand.1000\$y
	for (i in 1:(dim(df2.rand.1000)[1]-10)) {
		a = df2.rand.1000[df2.rand.1000\$x >= df2.rand.1000[i,]\$x & df2.rand.1000\$x <= df2.rand.1000[i+10-1,]\$x,]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2.rand.1000[i,]\$y2 = mean(a\$y)
			df2.rand.1000[i,]\$x2 = mean(a\$x)
		}
	}
	mins = seq(1,as.integer(df2.rand.1000[1,]\$x2)-1,10)
	maxs = seq(max(df2.rand.1000\$x2),dim(df.rand.1000)[2]-2,10)
	df2.rand.1000 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2.rand.1000)
	df2.rand.1000 = rbind(df2.rand.1000,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2.rand.1000 = data.frame(x=seq(1,dim(df.rand.1000)[2]), y=0, x2=seq(1,dim(df.rand.1000)[2]), y2=0);
}
 
# P2 % Conversion XY Plot
p2.rand.png.scale = p.png.scale
p2.rand.pdf.scale = p.pdf.scale

p2.rand.1000.png = 
	ggplot(df2.rand.1000,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.png.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.png.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

p2.rand.1000.pdf = 
	ggplot(df2.rand.1000,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.pdf.scale) +  # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.pdf.scale) + # SIZE 
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=5*p2.pdf.scale,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=5*p2.pdf.scale,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{Scale} = "

ratio1 = dim(df)[1] / (dim(df)[1] + 31.25 + 26.5625)
ratio2 = 31.25      / (dim(df)[1] + 31.25 + 26.5625)
ratio3 = 26.5625    / (dim(df)[1] + 31.25 + 26.5625)
mynrow = 2
totalheight = dim(df)[1] + 31.25
totalheight2 = dim(df)[1] + 31.25 + 26.5625
totalwidth  = dim(df)[2] * 0.125
totalratio  = c(ratio1, ratio2)
	if (file.exists(\"$curr_cluster_file\") & file.info(\"$curr_cluster_file\")\$size > 0) {
		totalheight = totalheight + 26.5625
		totalratio = c(totalratio, ratio3)
		mynrow = 3
	}

# Scaling

myscale = 4
totalheight  = totalheight * myscale
totalwidth = totalwidth * myscale
totalratio = totalratio * myscale
totalread_nopk = dim(df)[1] \%\% mywindow
totalheight_nopk = mywindow * myscale
totalheight_nopk_last = (totalread_nopk + 31.25) * myscale
totalratio_nopk_last  = c(totalread_nopk/(totalread_nopk+31.25+26.5625), 31.25 / (totalread_nopk+31.25+26.5625))
";

	my ($pngoutFolder, $pngoutFilename) = getFilename($pngout, "folderfull");
	my ($pdfoutFolder, $pdfoutFilename) = getFilename($pdfout, "folderfull");


	$R->{PNG} = "

# PNG
#png(type=\"cairo\",\"$pngout\",width=totalwidth,height=totalheight)
png(\"$pngout\",width=min(30000,totalwidth),height=min(30000,totalheight))
if (mynrow == 3) {
	grid.arrange(p.png,p2.png,p3.png,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p.png,p2.png,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.ALL.heatmap.png\"
#png(type=\"cairo\",pngout_heatmap_only,width=totalwidth,height=totalheight)
png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
print(p.heatmaponly)
dev.off()

# PNG all Conv
pngout_peak_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.ALL.c_conv.png\"
#png(type=\"cairo\",pngout_peak_all_c_conv,width=totalwidth,height=31.25*myscale)
png(pngout_peak_all_c_conv,width=min(30000,totalwidth),height=min(30000,31.25*myscale))
grid.arrange(p2.png)
dev.off()

";

	$R->{PNG_nopk} = "

# PNG
totalheight = (dim(df.rand.1000)[1] + 31.25) * myscale
png(\"$pngout\",width=min(30000,totalwidth),height=min(30000,totalheight))
#png(type=\"cairo\",\"$pngout\",width=totalwidth,height=min(30000,totalheight))
grid.arrange(p.rand.1000.png,p2.rand.1000.png,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df.rand.1000)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.ALL.heatmap.png\"
#png(type=\"cairo\",pngout_heatmap_only,width=min(30000,totalwidth),height=totalheight)
png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
print(p.heatmaponly.rand.1000)
dev.off()

# PNG all Conv
pngout_nopk_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.ALL.c_conv.png\"
png(pngout_nopk_all_c_conv,width=min(30000,totalwidth),height=31.25*myscale)
#png(type=\"cairo\",pngout_nopk_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2.png)
dev.off()

";

	$R->{PNG_nopk_rand_100} = "

pngout_nopk_rand_100 = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.png\"
# PNG
totalheight = (dim(df.rand.100)[1] + 31.25) * myscale
#png(type=\"cairo\",pngout_nopk_rand_100,width=totalwidth,height=totalheight)
png(pngout_nopk_rand_100,width=min(30000,totalwidth),height=min(30000,totalheight))
grid.arrange(p.rand.100.png,p2.rand.100.png,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df.rand.100)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.heatmap.png\"
#png(type=\"cairo\",pngout_heatmap_only,width=totalwidth,height=totalheight)
png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
print(p.heatmaponly.rand.100)
dev.off()

# PNG all Conv
pngout_nopk_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.c_conv.png\"
#png(type=\"cairo\",pngout_nopk_all_c_conv,width=totalwidth,height=31.25*myscale)
png(pngout_nopk_all_c_conv,width=min(30000,totalwidth),height=31.25*myscale)
grid.arrange(p2.rand.100.png)
dev.off()

";

	$R->{PNG_nopk_ALL} = "

for (i in seq(1,as.integer(dim(df)[1] / mywindow) + 1)) {
	pngout_nopk = \"$pngoutFolder/ALL/$pngoutFilename\"
	pngout_nopk = gsub(\"^(.+).png\$\",\"\\\\1\",pngout_nopk,perl=T)
	pngout_nopk = paste(pngout_nopk,\".\",i,\".png\",sep=\"\")
	if (i == max(seq( 1,as.integer(dim(df)[1]/mywindow) + 1 ))) {
		currtotalheight = totalheight_nopk_last
		currtotalratio = totalratio_nopk_last
		print(paste(i,currtotalheight,currtotalratio))
#		png(type=\"cairo\",pngout_nopk,width=totalwidth,height=currtotalheight)
		png(pngout_nopk,width=min(30000,totalwidth),height=min(30000,currtotalheight))
		grid.arrange(plot_list[[i]],p2,ncol=1,nrow=mynrow,heights=currtotalratio)
		dev.off()
	} else {
		currtotalheight = totalheight_nopk
		print(paste(i,currtotalheight))
	#	png(type=\"cairo\",pngout_nopk,width=totalwidth,height=currtotalheight)
		png(pngout_nopk,width=min(30000,totalwidth),height=min(30000,currtotalheight))
		grid.arrange(plot_list[[i]],ncol=1)
		dev.off()
	}
}


";

my $PDFSCALE = 200;
	$R->{PDF} = "

# PDF
currheight = totalheight / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdf(\"$pdfout\",width=min(30000,currwidth),height=min(30000,currheight))
if (mynrow == 3) {
	grid.arrange(p.pdf,p2.pdf,p3.pdf,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p.pdf,p2.pdf,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

# PDF HEATMAP ONLY
currheight = dim(df)[1] * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_heatmap_only = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.heatmap.pdf\"
pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
print(p.heatmaponly)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_peak_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_peak_all_c_conv,width=min(30000,currwidth),height=min(30000,currheight))
grid.arrange(p2.pdf)
dev.off()

";

	$R->{PDF_nopk} = "

# PDF
currheight = (dim(df.rand.1000)[1] + 31.25) * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdf(\"$pdfout\",width=min(30000,currwidth),height=min(30000,currheight))
grid.arrange(p.rand.1000.pdf,p2.rand.1000.pdf,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()


# PDF HEATMAP ONLY
currheight = dim(df.rand.1000)[1] * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_heatmap_only = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.heatmap.pdf\"
pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
print(p.heatmaponly.rand.1000)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=min(30000,currwidth),height=min(30000,currheight))
grid.arrange(p2.pdf)
dev.off()

";

	$R->{PDF_nopk_rand_100} = "
currheight = (dim(df.rand.100)[1] + 31.25) * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_rand_100 = \"$pdfoutFolder/ALL/$pdfoutFilename.RAND.100.pdf\"

# PDF
pdf(pdfout_nopk_rand_100,width=min(30000,currwidth),height=min(30000,currheight))
grid.arrange(p.rand.100.pdf,p2.rand.100.pdf,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PDF HEATMAP ONLY
currheight = dim(df.rand.100)[1] * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_heatmap_only = \"$pdfoutFolder/ALL/$pdfoutFilename.RAND.100.heatmap.pdf\"
pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
print(p.heatmaponly.rand.100)
dev.off()


# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.RAND.100.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=min(30000,currwidth),height=min(30000,currheight))
grid.arrange(p2.rand.100.pdf)
dev.off()

";

	$R->{PDF_nopk_ALL} = "

for (i in seq(1,as.integer(dim(df)[1] / mywindow) + 1)) {
	pdfout_nopk = \"$pdfoutFolder/ALL/$pdfoutFilename\"
	pdfout_nopk = gsub(\"^(.+).pdf\$\",\"\\\\1\",pdfout_nopk,perl=T)
	pdfout_nopk = paste(pdfout_nopk,\".\",i,\".pdf\",sep=\"\")
	if (i == max(seq( 1,as.integer(dim(df)[1]/mywindow) + 1 ))) {
		currtotalwidth  = totalwidth / $PDFSCALE
		currtotalheight = totalheight_nopk_last / $PDFSCALE
		currtotalratio = totalratio_nopk_last
		pdf(pdfout_nopk,width=currtotalwidth,height=min(30000,currtotalheight))
		grid.arrange(plot_list[[i]],p2,ncol=1,nrow=mynrow,heights=currtotalratio)
		dev.off()
	} else {
		currtotalwidth = totalwidth / $PDFSCALE
		currtotalheight = totalheight_nopk / $PDFSCALE
		pdf(pdfout_nopk,width=currtotalwidth,height=min(30000,currtotalheight))
		grid.arrange(plot_list[[i]],ncol=1)
		dev.off()
	}
}

";

	return $R;
}


__END__
