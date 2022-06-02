#!/usr/bin/perl
# footLoop pipeline Version 1.2
# Copyright (C) 2019 Stella Regina Hartono
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The license can be found at https://www.gnu.org/licenses/gpl-3.0.en.html. 
# By downloading or using this software, you agree to the terms and conditions of the license. 

use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars   qw($opt_v $opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f $opt_l $opt_e);
my @opts = qw($opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_l $opt_e);
getopts("vr:g:i:n:L:x:y:q:HhZFpl:e");

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
my $footLoopScriptsFolder = dirname(dirname abs_path $0);
$footLoopScriptsFolder .= "/footLoop";
my $cmd = "$footLoopScriptsFolder/check_software.pl | tail -n 12";
my @version = `$cmd`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

##################
# 0. Check Sanity #
##################

my %OPTS  = ('r' => $opt_r, 'g' => $opt_g, 'i' => $opt_i, 'n' => $opt_n, 
				 'L' => $opt_L, 'x' => $opt_x, 'y' => $opt_y, 'p' => $opt_p, 
				 'q' => $opt_q, 'F' => 'NONE', 'Z' => 'NONE',
				 'p' => 'NONE', 'q' => $opt_q, 'l' => $opt_l);
my %OPTS2 = ('p' => $opt_p, 'Z' => $opt_Z, 'F' => $opt_F);

check_sanity(\%OPTS, \%OPTS2);

###################
# 1. Define Input #
###################

my ($readFile, $genomeFile, $geneIndexFile, $outDir) = getFullpathAll($opt_r, $opt_g, $opt_i, $opt_n);
my ($readFilename)  = getFilename($opt_r, "full");
my ($geneIndexName) = getFilename($geneIndexFile);
my $minMapQ    = (not defined($opt_q)) ? 0  : $opt_q;
$opt_L = "95p" if not defined $opt_L;
my ($minReadL)   = $opt_L =~ /p$/i ? $opt_L =~ /^(.+)p$/i : $opt_L;
my $bufferL    = (not defined($opt_x)) ? 0  : $opt_x;
my $bufferR    = (not defined($opt_y)) ? 0  : $opt_y;
my $readName   = getFilename($readFile, "full");
my $bismarkOpt = (defined $opt_Z) ? "--bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8" : "--bowtie2 --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3";
my $samFile	   = ($readFile =~ /.f(ast)?q(.gz)?$/) ? $outDir .  "/$readName\_bismark_bt2.bam" : $readFile;
my $uuid       = getuuid();
my $date       = getDate();
my $origDir    = $outDir . "/.0_Orig";
my ($label)    = defined $opt_l ? $opt_l : $origDir =~ /PCB[\-_0]*\d+/i ? $origDir =~ /(PCB[\-_0]*\d+)/i : $outDir;
if (defined $opt_l) {
	$label = $opt_l;
}
elsif ($outDir =~ /PCB\d+/i) {
   ($label) = $outDir =~ /(PCB\d+)/i;
   $label = uc($label);
   $label =~ s/PCB0(\d+)/PCB$1/ if $label =~ /PCB0\d+/;
}
if (not defined $label and $outDir =~ /PCB[0\-_]*\d+/i) {
   ($label) = $outDir =~ /(PCB[0\-_]*\d+)/i;
   $label = uc($label);
   $label =~ s/PCB[0\-_]+(\d+)/PCB$1/g if $label =~ /PCB[0\-_]+\d+/;
}
if (not defined $label and $label =~ /PCB/) {
   ($label) = $outDir =~ /(PCB.{1,5})\/?/;
}
if (not defined $label) {
   die "Please make sure your output folder (-n) contain PCB(number) e.g. PCB12: -n 180202_PCB12_footloop_output\n\n";
}

if (not defined $opt_l and $outDir =~ /_BC\d+_PFC\d+_\w+\./i) {
   my ($label2) = $outDir =~ /_(BC\w+)\./i;
   $label = "$label\_$label2";
   $label = uc($label);
   $label =~ s/PCB[0\-_]+(\d+)/PCB$1/g;
}
	$OPTS{l} = $label;

if ($outDir =~ /PCB\d+_bcBC\d+_.+desc\w+/) {
	($label) = $outDir =~ /(PCB\d+_bcBC\d+_.+desc[a-zA-Z0-9]+)_?/;
}

print date() . "FATAL ERROR ON LABEL -l\n" and die if not defined $label;

# Make directory
makedir($outDir);

# make .LABEL
open (my $outLabel, ">", "$outDir/.LABEL") or die "Failed to write to $outDir/.LABEL: $!\n";
print $outLabel "$label\n";
close $outLabel;
my $STEP = 0;

# Make log file
my $logFile = "$outDir/logFile.txt";
open(my $outReadLog, ">", "$outDir/.PARAMS") or print "\n${LRD}FATAL!$N Failed to write to $LCY$outDir/.PARAMS$N: $!\n" and die;
open(my $outLog, '>', $logFile) or print "\n${LRD}FATAL!$N Failed to write to $LCY$logFile$N: $!\n" and die;
LOG($outReadLog, "footLoop.pl,uuid,$uuid\n");
LOG($outReadLog, "footLoop.pl,date,$date\n");
# Record all options

record_options(\%OPTS, \%OPTS2, $outReadLog, $version, $outLog);

#######################################
# 1. Preprocess Index and Fasta Files #
#######################################
($STEP) = LOGSTEP($outLog, "BEG", $STEP, 1, "Preprocessing gene index file and getting sequences\n");

# Add buffers to geneIndexFile ($footLoopScriptsFolder/lib/bedtools_bed_change.pl) and get their sequences using fastaFromBed, then get their info from fasta=$seqFile

my ($SEQ, $geneIndexHash, $seqFile, $bismark_folder) = parse_geneIndexFile($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog);

LOGSTEP($outLog);

###################
# 2. Runs Bismark #
###################
($STEP) = LOGSTEP($outLog, "BEG", $STEP, 1, "Creating bismark index and running bismark\n");#$N $bismarkOpt $bismark_folder $readFile\n");

# Make Bismark Index
make_bismark_index($seqFile, $bismark_folder, $bismarkOpt, $outLog);

# Run Bismark
($samFile) = run_bismark($readFile, $outDir, $samFile, $opt_F, $outReadLog, $outLog);

LOGSTEP($outLog);

################################
# 3. Fixing Sam File #
################################
($STEP) = LOGSTEP($outLog, "BEG", $STEP, 1, "Fix Sam File\n");#$N $bismarkOpt $bismark_folder $readFile\n");

# Do footLoop_2fixsam.pl
# - Determine strand of read based on # of conversion
# - Determine bad regions in read (indels) which wont be used in peak calling
($samFile, $origDir) = fix_samFile($samFile, $seqFile, $outReadLog, $outLog); #becomes .fixed
LOG($outReadLog, "footLoop.pl,origDir,$origDir\n");
my ($samFileName) = getFilename($samFile);

LOGSTEP($outLog);

################################
# 4. Parse and Filter Sam File #
################################
($STEP) = LOGSTEP($outLog, "BEG", $STEP, 1, "Parse and Filter Sam File\n");#$N $bismarkOpt $bismark_folder $readFile\n");

parse_samFile($samFile, $seqFile, $outReadLog, $outLog);

LOGSTEP($outLog);

#print_R_heatmap($SEQ);

###############
# SUBROUTINES #
###############

sub LOGSTEP {
	my ($outLog, $type, $STEP, $STEPCOUNT, $log2) = @_;
	if (not defined $type) {
		LOG($outLog, "$YW<==========$N\n\n");
		return;
	}
	($STEP) = LOG($outLog, "$YW==========>$N\nSTEP \$STEP\n", $STEP, $STEPCOUNT);
	if (defined $log2) {
		LOG($outLog, $log2);
	}
	return($STEP);
}

sub print_R_heatmap {
	my ($SEQ) = @_;
	

	my $box;
	foreach my $gene (sort keys %{$SEQ}) {
		my $Rscript = "$origDir/.$gene\_MakeHeatmap.R";
		open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
	
		# Breaks and color for R script
		# 0 is bad or not data (6)
		# 1 is non C/G
		# 4 is CH non conv
		# 5 is CG non conv
		# 6 is CH conv
		# 7 is CG conv
		# 8 is CH peak
		# 9 is CG peak
		
		my $breaks = "c(-1"	; my $colors = "c("               ; #BEGIN
		   $breaks .= ",0"	;    $colors .= "\"grey\""        ; # 0 is bad or no data (used to be 6)
		   $breaks .= ",1"	;    $colors .= ",\"white\""		 ; # 1 is non CG
		   $breaks .= ",4"	;    $colors .= ",\"cornsilk\""	 ; # 4 is CH non conv
			$breaks .= ",5"	;    $colors .= ",\"peachpuff\""  ; # 5 is CG non conv
			$breaks .= ",6"	;    $colors .= ",\"green4\""  	 ; # 6 is CH conv
			$breaks .= ",7"	;    $colors .= ",\"seagreen4\""  ; # 7 is CG conv
			$breaks .= ",8"	;    $colors .= ",\"red4\""    	 ; # 8 = PEAK Converted CH
			$breaks .= ",9"	;    $colors .= ",\"maroon4\""    ; # 9 = PEAK Converted CG
		   $breaks .= ")"    ;    $colors .= ")"               ; #END

		#Function to turn @seq into number and array for R
		my $seqR = "seq = c(";
		my @seq = @{$SEQ->{$gene}{seq}};
		for (my $s = 0; $s < @seq; $s++) {
			my $nuc = $seq[$s];
			#$nuc = $nuc =~ /A/i ? 11 : $nuc =~ /C/i ? 12 : $nuc =~ /G/i ? 13 : $nuc =~ /T/i ? 14 : 15;
			$nuc = $nuc =~ /A/i ? "seagreen4" : $nuc =~ /C/i ? "blue4" : $nuc =~ /G/i ? "saddlebrown" : $nuc =~ /T/i ? "red4" : "grey";
	
			$seqR .= "\"$nuc\"";
			$seqR .= $s == @seq - 1 ? ")" : ",";
		}
	
		my $labelz .= "p = p + geom_text(aes(x=min(dm\$x),label=id),hjust=0)";# if defined $opt_T;
		my $RBox = ""; my $RBoxPlot = "";
		if (defined $box) {
			if (defined $box->{$gene}) {
				my %RBox;
				foreach my $lines (sort @{$box->{$gene}}) {
					my ($chr, $beg, $end, $geneBox) = split("\t", $lines);
					push(@{$RBox{name}}, $geneBox);
					push(@{$RBox{beg}}, $beg);
					push(@{$RBox{end}}, $end);
					print "\t$geneBox, beg=$beg, end=$end\n";
				}
				if (defined $RBox{name}) {
					my $RBoxname = "RBoxname = c(\"" . join("\",\"", @{$RBox{name}}) . "\")";
					my $RBoxbeg = "RBoxbeg = c(" . join(",", @{$RBox{beg}}) . ")";
					my $RBoxend = "RBoxend = c(" . join(",", @{$RBox{end}}) . ")";
					$RBox = "
						mydm_xmin = min(dm\$x)
						mydm_xmax = max(dm\$x)
						$RBoxname
						$RBoxbeg
						$RBoxend
						if (length(RBoxbeg[RBoxbeg < mydm_xmin]) > 0) {
							RBoxbeg[RBoxbeg < mydm_xmin] = mydm_xmin
						}
						if (length(RBoxbeg[RBoxbeg > mydm_xmax]) > 0) {
							RBoxbeg[RBoxbeg > mydm_xmax] = mydm_xmax
						}
						if (length(RBoxend[RBoxend < mydm_xmin]) > 0) {
							RBoxend[RBoxend < mydm_xmin] = mydm_xmin
						}
						if (length(RBoxend[RBoxend > mydm_xmax]) > 0) {
							RBoxend[RBoxend > mydm_xmax] = mydm_xmax
						}
						RBox = data.frame(name=RBoxname,my_xmin=RBoxbeg,my_xmax=RBoxend, my_ymin=min(dm\$y),my_ymax = max(dm\$y))
					";
			
					$RBoxPlot = "+ geom_rect(data=RBox,aes(x=RBox\$my_xmin,y=RBox\$my_ymin,xmin=RBox\$my_xmin,xmax=RBox\$my_xmax,ymin=RBox\$my_ymin+0.5,ymax=RBox\$my_ymax-0.5),size=2,fill=NA,color=\"black\")";
				}
			}
		}
		print $out "
	
			.libPaths()
			library(\"GMD\")
			library(ggplot2)
			library(Cairo)
			#Sequence
			$seqR
	
	
			# print peak
			files = c()
			PDFout = c()
			PNGout = c()
			PEAKS = paste(files,\".PEAKS\",sep=\"\")
	
			info = file.info(files)
			files = rownames(info[info\$size != 0, ])
			for (i in 1:length(files)) 
			{
				if (file.exists(files[i])) 
				{
					myclust = TRUE
					if(length(grep(\"_NOPK\", files[i])) != 0) {
						myclust = FALSE
					}
					
					df = read.table(files[i],row.names=1,sep=\"\\t\")
					if (i > 2 & dim(df)[1] > 1000) {
						df = df[seq(1,1000),]
					}
	
					# in case there's only 1 peak:
					if (dim(df)[1] == 1) {
						df[2,] = df[1,]
					}
					if (length(df) > 0) {
	
	               rownames(df) = gsub(\"^SEQ_(\\\\w+)\$\",\"\\\\1\",rownames(df),perl=T)
	               ### Main ###
	               library(ggplot2)
	               library(reshape2)
	               library(grid)
	               h = hclust(d=dist(df))
	               df = df[h\$order,]
	               colnames(df) = seq(1,dim(df)[2])
	               df\$y = seq(1,dim(df)[1])
	               df\$id = rownames(df)
	               dm = melt(df,id.vars=c(\"id\",\"y\"))
	               colnames(dm) = c(\"id\",\"y\",\"x\",\"value\")
	               dm\$x = as.numeric(as.character(dm\$x))
	               dm\$y = as.numeric(as.character(dm\$y))
	               dm\$value = as.factor(as.character(dm\$value))
	
						if (file.exists(PEAKS[i])) {
	  	            	mypeaks = read.table(PEAKS[i]);
	   	            colnames(mypeaks) = c(\"id\",\"x\",\"xend\")
	      	         mydf = data.frame(id=rownames(df), y=df\$y)
	         	      mypeaks = merge(mypeaks, mydf, by=\"id\")
						} else {
							mypeaks=NULL
						}
					
	####### R BOX
	
					$RBox
	
	#############
	
	               p = ggplot(dm,aes(x,y)) +
	                   geom_tile(aes(fill=value)) + coord_fixed(ratio=4) +
	                   scale_fill_manual(values=c(\"0\"=\"cornsilk\",\"1\"=\"green4\",\"2\"=\"white\",\"3\"=\"peachpuff\",\"4\"=\"seagreen4\",\"5\"=\"maroon4\",\"6\"=\"grey\",\"9\"=\"red4\")) +
	                   scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	                   theme(line = element_blank(),
	                         axis.text = element_blank(),
	                         axis.title = element_blank()
	                   ) + ggtitle(paste(\"(\",dim(df)[1],\" peaks)\",sep=\"\")) $RBoxPlot
	################# R BOX PLOT
	
						if (length(mypeaks) > 0) {
							p = p + geom_rect(data=mypeaks,aes(xmin=mypeaks\$x,xmax=mypeaks\$xend,ymin=mypeaks\$y-0.5,ymax=mypeaks\$y+0.5),color=rgb(1,0,0),size=0.2,fill=NA); #+
									#geom_text(data=mypeaks,aes(x=mypeaks\$xend, y=mypeaks\$y, label=mypeaks\$id))
						}
						
						if (length(df) < 100) {
							$labelz
						}
	
	               gt <- ggplot_gtable(ggplot_build(p))
	               ge <- subset(gt\$layout, name == \"panel\")
	
	               ### PNG ###
	
	               print(paste(\"    \",i,\"B. Doing PNG of \",files[i],sep=\"\"))
	               png(type=\"cairo\",PNGout[i],height=4*2*dim(df)[1],width=2*dim(df)[2])
	               grid.draw(gt[ge\$t:ge\$b, ge\$l:ge\$r])
	               dev.off()
	
	
	
	               ### PNG ###
	               print(paste(i, \". Doing pdf of \",files[i],sep=\"\"))
	               pdf(PDFout[i])#,height=4*2*dim(df)[1],width=2*dim(df)[2]/4)
	               grid.draw(gt[ge\$t:ge\$b, ge\$l:ge\$r])
	               dev.off()
	               
					}
				}
			}
		";
		close $out;
		my $cmd = "R --no-save < $Rscript >> $logFile\n";
		print ($cmd);
		system($cmd) if -e $Rscript and -s $Rscript > 10;
		system("rm Rplots.pdf") if -e "Rplots.pdf";
		LOG($outLog, "\n${LPR}If there is not PDF made from step (8) then it's due to too low number of read/peak$N\n\n");
	}
}

sub fix_samFile {
	my ($samFile, $seqFile, $outReadLog, $outLog) = @_;
	LOG($outLog, "\n\ta. Fixing sam file $CY$samFile$N with $footLoopScriptsFolder/lib/footLoop_2fixsam.pl\n");
	my ($samMD5) = getMD5($samFile);
	my $origDir = "$outDir/.0_orig_$samMD5/";
	check_if_result_exist(["$origDir/.GOOD"], $outLog);
	makedir("$origDir") if not -d "$origDir";
	my $checkSam = 1;
	my ($samFileName) = getFilename($samFile, "full");
	my $samFileGZ = "$origDir/$samFileName.fixed.gz";
	$checkSam = 0 if not -e "$origDir/$samFileName.fixed" and not -e "$origDir/$samFileName.fixed.gz";
	makedir($origDir);
	if (defined $opt_F) {
		$checkSam = 0;
	}
	else {
		if (-e "$origDir/$samFileName.fixed.gz") {
			my ($samLineCount2) = linecount("$origDir/$samFileName.fixed.gz");
			my ($samLineCount1) = linecount($samFile);
			$checkSam = $samLineCount1 - 10 > $samLineCount2 ? 0 : 2;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: fixed sam file $LCY$origDir/$samFileName.fixed.gz$N exists but total row is less than total samFile $samFile row ($samLineCount1 - 500 > samFile.fixed.gz: $samLineCount2)!\n") if $checkSam == 0;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile::$LGN SUCCESS!!$N fixed sam file $LCY$origDir/$samFileName.fixed.gz$N exists (MD5=$LGN$samMD5$N) and total row $LGN($samLineCount2)$N >= total samFile row $LGN($samLineCount1 - 500)$N ($LCY$samFile$N)!\n") if $checkSam == 2;
		}
		if (-e "$origDir/$samFileName.fixed" and $checkSam == 0) {
			my ($samLineCount2) = linecount("$origDir/$samFileName.fixed");
			my ($samLineCount1) = linecount($samFile);
			$checkSam = $samLineCount1 - 10 > $samLineCount2 ? 0 : 1;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: .gz does not exist and fixed sam file $LCY$origDir/$samFileName.fixed$N exists but total row is less than total samFile $samFile row ($samLineCount1 - 500 > samFile.fixed: $samLineCount2)!\n") if $checkSam == 0;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile::$LGN SUCCESS!!$N fixed sam file $LCY$origDir/$samFileName.fixed$N exists (MD5=$LGN$samMD5$N) and total row $LGN($samLineCount2)$N >= total samFile row $LGN($samLineCount1 - 500)$N ($LCY$samFile$N)!\n") if $checkSam == 1;
		}
	}

	if ($checkSam == 0) {
		LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: fixed sam file $LCY$origDir/$samFileName.fixed$N or .gz does not exist!\n");
		LOG($outLog, "\t${YW}$footLoopScriptsFolder/lib/footLoop_2fixsam.pl -n $outDir -s $seqFile -o $origDir$N\n");
		system("$footLoopScriptsFolder/lib/footLoop_2fixsam.pl -n $outDir -o $origDir") == 0 or LOG($outLog, "Failed to run $footLoopScriptsFolder/lib/footLoop_2fixsam.pl -n $outDir -o $origDir: $!\n") and exit 1;
		LOG($outReadLog, "footLoop.pl,fix_samFile,$footLoopScriptsFolder/lib/footLoop_2fixsam.pl -n $outDir -o $origDir\n","NA");
		LOG($outLog, "\tgzip $origDir/$samFileName.fixed");
		system("gzip -f $origDir/$samFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $origDir/$samFileName.fixed: $!\n");
		$checkSam = 1;
	}
	else {
		LOG($outLog, "\t${LGN}WAS NOT RUN$N: ${YW}::: $footLoopScriptsFolder/lib/footLoop_2fixsam.pl -n $outDir -s $seqFile -o $origDir :::$N\n");
		# rm old (bad) .gz if it exists
		LOG($outLog, "\t/bin/rm $samFileGZ") if -e '$samFileGZ';
		system("/bin/rm $samFileGZ") == 0 or LOG($outLog, 'Failed to rm $samFileGZ: $!\n') if -e '$samFileGZ';
		# gzip the new .fixed
		LOG($outLog, "\tgzip $origDir/$samFileName.fixed");
		system("gzip -f $origDir/$samFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $origDir/$samFileName.fixed: $!\n");
	}

	# re-md5 samfile gz
	LOG($outLog, "\t${YW}$md5script $samFileGZ > $origDir/.$samFileName.fixed.gz.md5$N\n");
	if ($md5script eq "md5sum") {
		system("$md5script $samFileGZ > $origDir/.$samFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $samFileGZ > $origDir/.$samFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	if ($md5script eq "md5") {
		my ($md5res) = `$md5script $samFileGZ` =~ /^.+\= (\w+)$/; die "Failed to $md5script $samFileGZ: $!\n" if not defined $md5res;
		system("echo '$md5res $samFileGZ' > $origDir/.$samFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to $md5script $samFileGZ > $origDir/.$samFileName.fixed.gz.md5: $!\n") and exit 1;
	}
	LOG($outReadLog, "footLoop.pl,samFixedFile,$samFileGZ\n","NA");
	LOG($outReadLog, "footLoop.pl,samFixedFileMD5,$samMD5\n","NA");
	return($samFile, $origDir);
}

sub parse_samFile {
	my ($samFile, $seqFile, $outReadLog, $outLog) = @_;
	my ($samFileName) = getFilename($samFile);

	LOG($outLog, "\ta. Parsing sam file $CY$samFile$N and getting only high quality reads\n");
	open(my $notused, ">", "$outDir/.$readFilename.notused") or LOG($outLog, "Cannot open $outDir/.$readFilename.notused: $!\n") and exit 1;
	my $sam; 
	open($sam, $samFile) or LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!") and exit 1 if $samFile =~ /.sam$/;
	open($sam, "samtools view $samFile|") or LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!") and exit 1 if $samFile =~ /.bam$/;
	
	## Some stats
	my $linecount   = 0;
	my $samStats = makehash(['total','used','diffgene','lowq','badlength']);
	
	while(my $line = <$sam>) {
		$linecount ++;
		chomp($line);
		
		#####################
		# 1. Parse sam line #
		#####################
		my @arr = split("\t", $line);
		LOG($outLog, "\t$YW$samFile$N: Done $linecount\n") if $linecount % 5000 == 0;
	
		# a. Total column must be 14 or skipped (e.g. sam header)
		if (@arr < 14) {
			LOG($outLog, "$LGN$linecount$N: SAM header detected (#column < 14). LINE:\n\t$line\n");
			next;
		}
		
		my ($eval, $evalPrint) = myeval(\@arr);
		my ($read, $readStrand, $gene, $readPos, $mapQ, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $tags, $junk4, $readMethCall) = @arr;
		my $seqLen = length($seqs);
		$gene = uc($gene);
		check_sam_field($samFile, $outLog, @arr);
		LOG($outLog, "\t$samFile: length of gene $gene undef\n") and exit 1 if not defined($SEQ->{$gene}{geneL});

		$samStats->{total} ++;
		$samStats->{gene}{$gene} ++;
	
		# Filter out reads which gene isn't in genomeFile and indexFile
		if (not defined $SEQ->{$gene}) {
			LOG($outLog, "$LGN$linecount$N: read $LCY$read$N is mapped to gene $LGN$gene$N is not in geneIndex! LINE:\n\t$line\n") if not defined $SEQ->{$gene};
			next;
		}
		
		$SEQ->{$gene}{total} ++;
		my ($readname) = length($read) >= 20 ? $read =~ /(.{20})$/ : $read; $readname = "..." . $readname  if length($read) >= 20; $readname = $read if not defined $readname;
		LOG($outLog, "\tExample at read $samStats->{total}: name=$CY$readname$N\tstrand=$CY$readStrand$N\tchr/gene=$CY$gene$N\tpos=$CY$readPos$N\tmapQ=$CY$mapQ$N\n\n") if $samStats->{total} == 1;
		LOG($outLog, "\tDone $GN$samStats->{total}$N\n") if $samStats->{total} % 2000 == 0;

		# a. Filter out reads with mapping quality=$mapQ less than threshold quality=$opt_q
		if($mapQ < $minMapQ)	{
			($SEQ, $samStats) = filter_sam_read($SEQ, $samStats, $readname, $gene, 'lowq');
			print $notused "\t$CY$read$N quality ($CY$mapQ$N) is less than $CY$opt_q$N\n";
			$samStats->{lowq} ++;
			$SEQ->{$gene}{lowq} ++;
			next;
		}
		
		# b. Filter out reads with read length=$seqLen
		elsif (parse_cigar($cigar, "len") < $SEQ->{$gene}{minReadL}) {
			($SEQ, $samStats) = filter_sam_read($SEQ, $samStats, $readname, $gene, 'badlength');
			my $cigarLen = parse_cigar($cigar, "len");
			print $notused "\t$CY$read$N length of seq (cigarLen = $LGN$cigarLen$N, seqLen = $CY$seqLen$N) is less than $SEQ->{$gene}{minReadL} bp (length of original sequence is ($CY" . $SEQ->{$gene}{geneL} . "$N)!\n";
		}

		# c. If read strand is 0 or 16 store it
		elsif ($readStrand == 0 || $readStrand == 16)	{
			$SEQ->{$gene}{read}{$read} = $readStrand;
		}
		
		# d. Otherwise put it into notused
		else {
			print $notused "\t$CY$read$N isn't used for some reason\n";
		}
	}
	LOG($outLog, "\n\tb.Logging sam file $CY$samFile$N\n");
	log_samFile($SEQ, $samStats, $origDir, $outReadLog, $outLog);
}

sub check_sam_field {
	my ($samFile, $outLog, @arr) = @_;
	my $check = 0;
	for (my $i = 0; $i < @arr; $i++) {
		if ( not defined $arr[$i]) {
			LOG($outLog, "\t$samFile: column $i undefined\n");
			$check = 1;
		}
	}
	exit 1 if $check != 0;
}

sub log_samFile {
	my ($SEQ, $samStats, $origDir, $outReadLog, $outLog) = @_;
	foreach my $genez (sort keys %{$SEQ}) {
		my $outTXTFilePos  = "$origDir/$genez\_Pos.orig";
		my $outTXTFileNeg  = "$origDir/$genez\_Neg.orig";
		my $outTXTFileUnk  = "$origDir/$genez\_Unk.orig";
		open ($SEQ->{$genez}{outTXTPos}, ">", $outTXTFilePos);
		open ($SEQ->{$genez}{outTXTNeg}, ">", $outTXTFileNeg);
		open ($SEQ->{$genez}{outTXTUnk}, ">", $outTXTFileUnk);
	}
	
	my $skipped = 0; my ($passedFilterP, $passedFilterN) = (0,0);
	open (my $inSamFix, "zcat < $origDir/$samFileName.fixed.gz|") or LOG($outLog, "Failed to open $origDir/$samFileName.fixed.gz: $!\n") and exit 1;
	while (my $line = <$inSamFix>) {
		chomp($line);
		my ($read, $type, $oldStrand, $strand, $genez, $pos, $info) = split("\t", $line);
		if (not defined $SEQ->{$genez} or not defined $info) {
			LOG($outLog, "\tERROR in $LCY$origDir/$samFileName.fixed$N: gene=$genez but \$SEQ->{\$genez} is not defined!line=\n$line\n\n") if not defined $SEQ->{$genez};
			DIELOG($outLog, "\tERROR in $LCY$origDir/$samFileName.fixed$N: gene=$genez and seq genez is $SEQ->{$genez} but info is not defined!line=\n$line\n\n") if defined $SEQ->{$genez} and not defined $info;
			$skipped ++;
			next;
		}
		if (not defined $SEQ->{$genez}{read}{$read}) {
			next;
		}
		my ($CT0, $CC0, $GA0, $GG0, $CT1, $CC1, $GA1, $GG1) = split(",", $info);
	#	my $oldStrand = $SEQ->{$genez}{read}{$read};
		if ($type eq "6_BOTH" or $type eq "3_NONE") {
			print {$SEQ->{$genez}{outTXTUnk}} "$read\tBP\t$pos\n" if $strand eq 0;
			print {$SEQ->{$genez}{outTXTUnk}} "$read\tBN\t$pos\n" if $strand eq 16;
			$SEQ->{$genez}{unkpos} ++ if $strand == 0;
			$SEQ->{$genez}{unkneg} ++ if $strand == 16;
			$passedFilterP ++ if $strand == 0;
			$passedFilterN ++ if $strand == 16;
		}
		else {
			print {$SEQ->{$genez}{outTXTNeg}} "$read\tFN\t$pos\n" if $strand eq 16;
			print {$SEQ->{$genez}{outTXTPos}} "$read\tFP\t$pos\n" if $strand eq 0;
			$SEQ->{$genez}{pos} ++ if $strand == 0;
			$SEQ->{$genez}{neg} ++ if $strand == 16;
			$passedFilterP ++ if $strand == 0;
			$passedFilterN ++ if $strand == 16;
		}
		$samStats->{used} ++;
		$SEQ->{$genez}{used} ++;
		$SEQ->{$genez}{posneg} ++ if $strand eq 16 and $oldStrand eq 0;
		$SEQ->{$genez}{negpos} ++ if $strand eq 0 and $oldStrand eq 16;
	}
	close $inSamFix;
	LOG($outLog, "\t$LCY$origDir/$samFileName.fixed$N: skipped = $LGN$skipped$N\n");

	LOG($outReadLog, "footLoop.pl,read_passed_filter,header\ttotal\tpositive\tnegative\n");
	LOG($outReadLog, "footLoop.pl,read_passed_filter,record\t$samStats->{total}\t$passedFilterP\t$passedFilterN\n");

	print $outLog "
Reads that passed filters:
Positive: $passedFilterP
Negative: $passedFilterN
Total   : $samStats->{total};
	
Per Gene:
	";
	
	my $zero = "";
	foreach my $gene (sort keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		foreach my $key (@key) {
			$SEQ->{$gene}{$key} = 0 if not defined $SEQ->{$gene}{$key};
		}
	}
	foreach my $gene (sort {$SEQ->{$b}{total} <=> $SEQ->{$a}{total}} keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		my $outTXTFilePos  = "$origDir/$gene\_Pos.orig"; system("/bin/rm $outTXTFilePos") if -e $outTXTFilePos and -s $outTXTFilePos == 0;
		my $outTXTFileNeg  = "$origDir/$gene\_Neg.orig"; system("/bin/rm $outTXTFileNeg") if -e $outTXTFileNeg and -s $outTXTFileNeg == 0;
		my $outTXTFileUnk  = "$origDir/$gene\_Unk.orig"; system("/bin/rm $outTXTFileUnk") if -e $outTXTFileUnk and -s $outTXTFileUnk == 0;
	
		$zero .= "$gene ($SEQ->{$gene}{total})\n" and next if $SEQ->{$gene}{total} <= 10;
		my $gene2 = $SEQ->{$gene}{orig};
		my $text = "
- $gene (original name = $gene2):
Positive    = $SEQ->{$gene}{pos} (neg->pos = $SEQ->{$gene}{negpos})
Negative    = $SEQ->{$gene}{neg} (pos->neg = $SEQ->{$gene}{posneg})
UnkPos      = $SEQ->{$gene}{unkpos}
UnkNeg      = $SEQ->{$gene}{unkneg}
Used        = $SEQ->{$gene}{used}
Total       = $SEQ->{$gene}{total}
Too Short   = $SEQ->{$gene}{badlength}
Low Quality = $SEQ->{$gene}{lowq}
";
	LOG($outLog, $text);
	LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,header\tgene\ttotal\tused\tpositive\tnegative\tunkpos\tunkneg\ttooshort\tlowqual\n");
	LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,record\t$gene\t$SEQ->{$gene}{total}\t$SEQ->{$gene}{used}\t$SEQ->{$gene}{pos}\t$SEQ->{$gene}{neg}\t$SEQ->{$gene}{unkpos}\t$SEQ->{$gene}{unkneg}\t$SEQ->{$gene}{badlength}\t$SEQ->{$gene}{lowq}\n");
	LOG($outReadLog, "footLoop.pl,gene_skipped_low_read,$zero\n");
	}
	
	$zero = $zero eq "" ? "(None)\n" : "\n$zero\n";

	LOG($outLog, "
	(Genes that have <= 10 reads:) $LGN$zero$N
	${GN}SUCCESS$N: Total=$samStats->{total}, used=$samStats->{used}, Low Map Quality=$samStats->{lowq}, Too short=$samStats->{badlength}
	
	Output: ${LGN}$origDir$N\n
	");

}
sub filter_sam_read {
	my ($seq, $count, $readname, $gene, $reason) = @_;
	$count->{$reason} ++;
	$seq->{$gene}{$reason} ++;
	return($seq, $count);
}

sub get_lotsOfC {
	my ($seq) = @_;
	my @seq = $seq =~ /ARRAY/ ? @{$seq} : split("", $seq);
	my $len = @seq;
	my $data;
	my %nuc;
	foreach my $nuc (@seq[0..@seq-1]) {
		$nuc = uc($nuc);
		$nuc{$nuc} ++;
	}
	foreach my $nuc (sort keys %nuc) {
		$nuc = uc($nuc);
		next if $nuc{$nuc} < 10;
		next if $nuc !~ /^(C|G)$/i;
		while ($seq =~ /${nuc}{6,$len}/ig) {
			my ($prev, $curr, $next) = ($`, $&, $');
			my ($beg, $end) = (length($prev), length($prev) + length($curr));
			$data .= ",$nuc;$beg;$end";
			my ($prev0) = length($prev) > 2 ? $prev =~ /(..)$/ : "NOPREV";
			my ($next0) = length($next) > 2 ? $next =~ /^(..)/ : "NONEXT";
		}
	}
	$data =~ s/^,// if defined $data;
	return($data);
}

sub make_bismark_index {
	my ($geneIndexFa, $bismark_folder, $bismarkOpt, $outLog) = @_;
	LOG($outLog, "\n\ta. Running bismark_genome_preparation$N --bowtie2 $bismark_folder$N\n");
	my $run_boolean = "\t${LGN}WAS NOT RUN$N:${YW} ";
	my $cmd = "bismark_genome_preparation --bowtie2 $bismark_folder > $bismark_folder/LOG.txt 2>&1";
	if ($md5script eq "md5sum") {
		$cmd .= " && $md5script $geneIndexFa > $bismark_folder/Bisulfite_Genome/.$md5script";
	}
	my ($check, $md5sum, $md5sum2) = (0);
	my $bismark_folder_exist = (-d "$bismark_folder/Bisulfite_Genome/" and -e "$bismark_folder/Bisulfite_Genome/.MD5SUM") ? 1 : 0;
	if ($bismark_folder_exist == 1) {
		LOG($outLog, "\tOlder bismark folder $CY$bismark_folder/BIsulfite_Genome/$N exist! Checking MD5 if they're the same as current fasta file.\n");
		($md5sum)  = `cat $bismark_folder/Bisulfite_Genome/.MD5SUM` =~ /^(\w+)[\t ]/;
		($md5sum2) = getMD5($geneIndexFa);
		$bismark_folder_exist = 0 if $md5sum ne $md5sum2;
	}
	$run_boolean = "" if $bismark_folder_exist == 0;
	if ($bismark_folder_exist == 0) {
		LOG($outLog, "\tEither bismark folder didn't exist or older bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n") if defined $md5sum;
		system($cmd) == 0 or die "Failed to run bismark genome preparation: $!\n";
	   if ($md5script eq "md5") {
	      my ($md5res) = `$md5script $geneIndexFa` =~ /^.+\= (\w+)$/; die "Failed to $md5script $geneIndexFa: $!\n" if not defined $md5res;
	      system("echo '$md5res $geneIndexFa' > $bismark_folder/Bisulfite_Genome/.$md5script") == 0 or LOG($outLog, "Failed to $md5script $geneIndexFa: $!\n");
	   }
	}
	else { 
		LOG($outLog, "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n");
	}
	LOG($outLog, "${run_boolean} $YW ::: $cmd :::$N\n");
}

sub run_bismark {
	my ($readFile, $outDir, $mybam, $opt_F, $outReadLog, $outLog) = @_;
	my $MAP = "";
	LOG($outLog, "\n\tb. Running bismark\n");
	$mybam =~ s/.fq.gz_bismark_bt2/_bismark_bt2/;
	my $bamFile = $mybam; 
	my $ext = "bam";
	
	my ($mybamFilename) = getFilename($mybam, "full");
	
	my $run_boolean = "\n\t${LGN}WAS NOT RUN$N:${YW} ";

	if (-e $mybam and not -e "$outDir/$mybamFilename") {
		system("/bin/ln -s $mybam $outDir/$mybamFilename") == 0 or LOG($outLog, "Failed to /bin/ln $mybam $outDir/$mybamFilename: $!\n") and exit 1;
	}
	if (defined $opt_F or not -e $mybam) {
		$run_boolean = "\n$YW\t ";
		if ($opt_p) {
			my $outFolder = $outDir . "/.bismark_paralel/";
			if (-d $outFolder) {
				my @files = <$outFolder/*>;
				if (@files != 0) {
					system("/bin/rm $outFolder/*");
				}
			}
			makedir($outFolder) if not -d $outFolder;
			LOG($outLog, "\t  Splitting $CY$readFile$N by 1000 sequences!\n\n####### SPLITRESULT LOG ########");	
			my $splitresult = `SplitFastq.pl -i $readFile -o $outFolder -n 1000`;
			LOG($outReadLog, "footLoop.pl,run_bismark,SplitFastq.pl -i $readFile -o $outFolder -n 1000\n");

			LOG($outLog, "$splitresult\n");	
			LOG($outLog, "####### SPLITRESULT LOG #######\n\n\t  Running bismark in paralel!\n");
			#my $result = system("run_script_in_paralel2.pl -v \"srun -p high --mem 8000 bismark -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\" $outFolder .part 20");
			my $result = system("run_script_in_paralel2.pl -v \"srun -p high --mem 8000 bismark --non_directional -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\" $outFolder .part 20");
			#LOG($outReadLog, "footLoop.pl,run_bismark,\"run_script_in_paralel2.pl -v \\\"srun -p high --mem 8000 bismark -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\\\" $outFolder .part 20");
			LOG($outReadLog, "footLoop.pl,run_bismark,\"run_script_in_paralel2.pl -v \\\"srun -p high --mem 8000 bismark --non_directional -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\\\" $outFolder .part 20");
			my @partSam = <$outFolder/*.part_bismark*.$ext>; my $totalPartSam = @partSam;
			LOG($outLog, "\t  All part.$ext has been made (total = $totalPartSam). Now making $CY$mybam$N and then removing the part sam\n");
			my @HEADER; my @REPORT;
			for (my $p = 0; $p < @partSam; $p++) {
				my $partSam = $partSam[$p];
				print "\t\tPutting $partSam into $mybam and removing it!\n";
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mybam") == 0 or die "Failed to cat $partSam: $!\n" if $p == 0;
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mybam") == 0 or die "Failed to cat $partSam: $!\n" if $p != 0;
				LOG($outReadLog, "footLoop.pl,run_bismark,cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mybam") if $p == 0;
				LOG($outReadLog, "footLoop.pl,run_bismark,cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mybam") if $p != 0;
				my ($bismark_report) = $partSam . "_SE_report.txt";
				if (not -e $bismark_report) {
					($bismark_report) =~ /^(.+).$ext/;
					$bismark_report .= "_SE_report.txt";
				}

				my ($header, $report) = parse_bismark_report($bismark_report, $outLog);
				my @header = @{$header};
				my @report = @{$report};
				@HEADER = @header if $p == 0;
				for (my $q = 0; $q < @header; $q++) {
					die "undefined $q header\n" if not defined $header[$q];
					die if $header[$q] ne $HEADER[$q];
					$REPORT[$q] += $report[$q] if $header[$q] !~ /^perc_/;
					$REPORT[$q] += $report[0]*$report[$q] if $header[$q] =~ /^perc_/;
				}
			}
			for (my $q = 0; $q < @HEADER; $q++) {
				$REPORT[$q] = int(100*$REPORT[$q]/$REPORT[0]+0.5)/100 if $HEADER[$q] =~ /^perc_/;
			}
			$MAP  = "footLoop.pl,map," . "header\tlabel\t" . join("\t", @HEADER) . "\tfootLoop_outDir\tuuid\n";
		   $MAP .= "footLoop.pl,map," . "record\t$label\t" . join("\t", @REPORT) . "\t$outDir\t$uuid\n";
		}
		else {
			LOG($outLog, "\t  bismark --non_directional -o <outDir> $LCY$bismarkOpt$N <bismark_folder> <readFile> > <outDir/.bismark_log> 2>&1\n");
			#LOG($outLog, "\t  bismark -o <outDir> $LCY$bismarkOpt$N <bismark_folder> <readFile> > <outDir/.bismark_log> 2>&1\n");
			LOG($outReadLog, "footLoop.pl,bismark,bismark --non_directional -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1\n","NA");
			#LOG($outReadLog, "footLoop.pl,bismark,bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1\n","NA");
			my $result = system("bismark --non_directional -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");
			#my $result = system("bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");

			if ($result != 0) {
				LOG($outLog, "\t\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_folder\n");
				make_bismark_index($seqFile, $bismark_folder, $bismarkOpt, $outLog);
				system("bismark_genome_preparation --bowtie2 $bismark_folder") == 0 or die "Failed to run bismark genome preparation: $!\n";
				LOG($outReadLog, "footLoop.pl,bismark,bismark_genome_preparation --bowtie2 $bismark_folder","NA");
				system("bismark --non_directional -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
				LOG($outReadLog, "footLoop.pl,bismark,bismark --non_directional -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");
				#LOG($outReadLog, "footLoop.pl,bismark,bismark_genome_preparation --bowtie2 $bismark_folder","NA");
				#system("bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
				#LOG($outReadLog, "footLoop.pl,bismark,bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");
			}
			LOG($outLog, "\t${GN}SUCCESS$N: Output $mybam\n");
			my ($bismark_report) = $mybam =~ /^(.+).$ext/; $bismark_report .= "_SE_report.txt";
			my ($header, $report, $MAPTEMP) = parse_bismark_report($bismark_report, $outLog);
			$MAP = $MAPTEMP;
		}
	}
	else {
		LOG($outLog, "\t${GN}SUCCESS$N: Output already exist: $CY$mybam$N\n");
		my ($bismark_report) = "$mybam" =~ /^(.+).$ext/; $bismark_report .= "_SE_report.txt";
		my ($header, $report, $MAPTEMP) = parse_bismark_report($bismark_report, $outLog);
		$MAP = $MAPTEMP;
	}
	LOG($outLog, "${run_boolean}::: bismark $bismarkOpt $bismark_folder $readFile :::$N\n");

	# print to $outReadLog
	LOG($outReadLog, "footLoop.pl,samFile,$outDir/$mybamFilename\n","NA");
	LOG($outReadLog, $MAP,"NA");

	return("$outDir/$mybamFilename");
}

sub parse_bismark_report {
	my ($bismark_report_file, $outLog) = @_;
	open (my $inz, "<", $bismark_report_file) or DIELOG($outLog, "Failed to read from $LCY$bismark_report_file$N: $!\n");
	my @report; my @header;
	while (my $line = <$inz>) {
		chomp($line);
		if ($line =~ /^Sequences analysed in total:/) {
			my ($num) = $line =~ /in total:[ \t]+(\d+)[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get total read from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "total_read");
		}
		if ($line =~ /^Mapping efficiency/) {
			my ($num) = $line =~ /efficiency:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Percent Mapping Efficiency from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "perc_mapped_read");
		}
		if ($line =~ /^Sequences which were discarded because genomic sequence could not be extracted/) {
			my ($num) = $line =~ /extracted:[ \t]+(\d+\.?\d*)[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Cannot Extract Chr from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "cannot_extract_chr");
		}
		if ($line =~ /^Total number of C's analysed:/) {
			my ($num) = $line =~ /Total number of C's analysed:[ \t]+(\d+\.?\d*)[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total C analysed from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "total_C");
		}
		if ($line =~ /^C methylated in CpG context:/) {
			my ($num) = $line =~ /C methylated in CpG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CpG from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "perc_meth_CpG");
		}
		if ($line =~ /^C methylated in CHG context:/) {
			my ($num) = $line =~ /C methylated in CHG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHG from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "perc_meth_CHG");
		}
		if ($line =~ /^C methylated in CHH context:/) {
			my ($num) = $line =~ /C methylated in CHH context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHH from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "perc_meth_CHH");
		}
		if ($line =~ /^C methylated in Unknown context:/) {
			my ($num) = $line =~ /C methylated in Unknown context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in Unknown from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "perc_meth_Unknown");
		}
	}
	close $inz;
	my $MAP  = "footLoop.pl,map," . "header\tlabel\t" . join("\t", @header) . "\tfootLoop_outDir\tuuid\n";
	   $MAP .= "footLoop.pl,map," . "record\t$label\t" . join("\t", @report) . "\t$outDir\t$uuid\n";
	return(\@header, \@report, $MAP);

}


sub get_geneIndex_fasta {
	my ($geneIndexFile, $outDir, $logFile, $outLog) = @_;
	my $geneIndexName = getFilename($geneIndexFile);
	my $run_boolean = "\t${LGN}WAS NOT RUN$N:${YW} ";
	my $geneIndexFileNew = "$outDir/$geneIndexName\_$bufferL\_$bufferR\_bp.bed";
	LOG($outLog, "\n\ta. Transforming $LCY $geneIndexFile$N into$LCY $geneIndexFileNew$N\n");
	if ($bufferL eq 0 and $bufferR eq 0) {
		LOG($outLog, "$YW\t::: /bin/cp $geneIndexFile $geneIndexFileNew :::$N\n");
		system("/bin/cp $geneIndexFile $geneIndexFileNew") == 0 or LOG($outLog, "\tfootLoop.pl::get_geneIndex_fasta: Failed to$YW /bin/cp $geneIndexFile $geneIndexFileNew$N: $!\n") and exit 1;
		LOG($outLog, "\t${GN}SUCCESS$N: Created index file $LGN$geneIndexFileNew$N from $LCY$geneIndexFile$N\n");
	}
	else {
		LOG($outLog, "$YW\t::: $footLoopScriptsFolder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1 :::$N\n") == 0 or LOG($outLog, "Failed to get (beg $bufferL, end $bufferR) bp of $geneIndexFile!\n") and exit 1;
		system("$footLoopScriptsFolder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1") == 0 or LOG($outLog, "\tfootLoop.pl::get_geneIndex_fasta: Failed to$YW $footLoopScriptsFolder/lib/bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1$N\n: $!\n") and exit 1;
		LOG($outLog, "\t${GN}SUCCESS$N: Created index file $LGN$geneIndexFileNew$N from bedtools bed change of $LCY$geneIndexFile$N\n");
	}
	return($geneIndexFileNew);
}

sub parse_geneIndexFile {
	my ($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog) = @_;

	my $geneIndexlinecount = 0;
	open (my $geneIndexInCheck, "<", $geneIndexFile) or DIELOG($outLog, "\nCannot read from geneIndexFile $geneIndexFile: $!\n");
	while (my $line = <$geneIndexInCheck>) {
		chomp($line); $geneIndexlinecount ++;
		next if $line =~ /^#/;
		my @arr = split("\t", $line);
		if (@arr < 6) {
			DIELOG($outLog, "\nERROR: geneIndexFile (-i $geneIndexFile) has to be a 6 column bed format! Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[1] !~ /^\d+$/) {
			DIELOG($outLog, "\nERROR: geneIndexFile (-i $geneIndexFile) column 1 isn't integer! ($arr[1]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[2] !~ /^\d+$/) {
			DIELOG($outLog, "\nERROR: geneIndexFile (-i $geneIndexFile) column 2 isn't integer! ($arr[2]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[4] !~ /^\-?\d+\.?\d*e?\-?\d*\.?\d*$/i) {
			DIELOG($outLog, "\nERROR: geneIndexFile (-i $geneIndexFile) column 4 isn't numeric ($arr[4]) Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
		if ($arr[5] !~ /^[\+\-]$/) {
			DIELOG($outLog, "\nERROR: geneIndexFile (-i $geneIndexFile) column 6 isn't strand (- or +)! Offending line linecount = $geneIndexlinecount\n\n$line\n\n");
		}
	}
	close $geneIndexInCheck;

	$geneIndexFile = get_geneIndex_fasta($geneIndexFile, $outDir, $logFile, $outLog);
	LOG($outReadLog, "footLoop.pl,geneIndexFile,$geneIndexFile\n","NA");
	my $geneIndex; my $linecount = 0;
	open (my $geneIndexIn, "<", $geneIndexFile) or die "Cannot read from $geneIndexFile: $!\n";
	LOG($outLog, "\t\t${GR}From geneIndexFile:$N\n");
	while (my $line = <$geneIndexIn>) {
		chomp($line); $linecount ++;
		my ($chr, $beg, $end, $gene) = split("\t", $line);
		$gene = uc($gene);
		LOG($outLog, "\t\t$GR$linecount: gene=$gene,beg=$beg,end=$end,length=" . ($end-$beg) . "$N\n");
		$geneIndex->{$gene} = $beg;
	}
	close $geneIndexIn;
	LOG($outLog, "\t${GN}SUCCESS$N: Parsed gene coordinates from index file $LGN$geneIndexFile$N\n");

	my $geneIndexName = getFilename($geneIndexFile, "full");
	my $geneIndexFaTemp = $outDir . "/.geneIndex/$geneIndexName.fa";
	makedir($geneIndexFaTemp, 1);

	LOG($outLog, "\n\tb. Getting fasta sequence from bed file$LCY $geneIndexFile$N genome file$LCY $genomeFile$N\n");

	my $cmd = "fastaFromBed -fi $genomeFile -bed $geneIndexFile -fo $geneIndexFaTemp -name";
	LOG($outLog, "$YW\t::: $cmd :::$N\n");
	system("$cmd") == 0 ? LOG($outLog, "\t${GN}SUCCESS$N: Output: geneIndexFaTemp=$CY$geneIndexFaTemp$N\n") : die "Failed to run bedtools: $!\n";

	$geneIndexFaTemp = uppercaseFasta($geneIndexFaTemp, $outLog);
	my ($geneIndexFaMD5, $temp, $geneIndexFaTempMD5File)  = getMD5($geneIndexFaTemp);
	LOG($outLog, "\tFatal error as MD5 of $geneIndexFaTemp is NA..(GENE=$geneIndexFaTemp, MD5=$geneIndexFaMD5)\n") and exit 1 if $geneIndexFaMD5 eq "NA";
	my $geneIndexFa        = $outDir . "/.geneIndex/$geneIndexFaMD5/$geneIndexName.fa";
	my $geneIndexFaMD5File = $outDir . "/.geneIndex/$geneIndexFaMD5/.$geneIndexName.fa.md5";
	if (not -d $geneIndexFaMD5 or (-d $geneIndexFaMD5 and not -e $geneIndexFa)) {
		makedir("$outDir/.geneIndex/$geneIndexFaMD5");
		system("/bin/mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/") == 0 or LOG($outLog, "Failed to mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/: $!\n") and exit 1;
		system("/bin/mv $geneIndexFaTempMD5File $outDir/.geneIndex/$geneIndexFaMD5/") == 0 or LOG($outLog, "Failed to mv $geneIndexFaTemp $outDir/.geneIndex/$geneIndexFaMD5/: $!\n") and exit 1;
	}
	LOG($outLog, "\nERROR: geneIndexHash is not defined!\n") and exit 1 if not defined $geneIndex;
	LOG($outLog, "\nERROR: geneIndexFa is not defined!\n")   and exit 1 if not defined $geneIndexFa;
	my $SEQ = parse_fasta($geneIndexFa, $outLog, $minReadL);
	LOG($outReadLog, "footLoop.pl,-g,$genomeFile\n","NA");
	LOG($outReadLog, "footLoop.pl,seqFile,$geneIndexFa\n","NA");
	return ($SEQ, $geneIndex, $geneIndexFa, "$outDir/.geneIndex/$geneIndexFaMD5/");
	#my ($SEQ, $geneIndexHash, $seqFile, $bismark_folder) = parse_geneIndexFile($geneIndexFile, $outDir, $outLog, $seqFile, $minReadL);
}

sub parse_fasta {
	my ($seqFile, $outLog, $minReadL) = @_;
	LOG($outLog, "\n\tc. Parsing in gene sequence and infos from seqFile=$CY$seqFile$N\n");
	open(my $SEQIN, "<", $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
	my $fasta = new FAlite($SEQIN);
	my $linecount = 0;
	LOG($outLog, "\t\t${GR}From fasta file:$N\n");
	while (my $entry = $fasta->nextEntry()) {
		$linecount ++;
		my $def = $entry->def; $def =~ s/^>//;
		my $gene = uc($def);
		my $seqz = uc($entry->seq);
	   $SEQ->{$gene}{seq}       = [split("", $seqz)];
		$SEQ->{$gene}{geneL}     = scalar(@{$SEQ->{$gene}{seq}});
		$SEQ->{$gene}{minReadL}  = (defined $minReadL and $opt_L =~ /p$/i) ? int(0.5+$SEQ->{$gene}{geneL} * $minReadL / 100) : $minReadL;
		$SEQ->{$gene}{total}     = 0;
		$SEQ->{$gene}{badlength} = 0;
		$SEQ->{$gene}{lowq}      = 0;
		$SEQ->{$gene}{used}      = 0;
		$SEQ->{$gene}{pos}       = 0;
		$SEQ->{$gene}{neg}       = 0;
		$SEQ->{$gene}{orig}      = $def;
		LOG($outLog, "\t\t$GR$linecount:gene=$gene,length=$SEQ->{$gene}{geneL}$N\n");
	}
	close $SEQIN;
	LOG($outLog, "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n");
	return($SEQ);
}

sub uppercaseFasta {
	my ($fastaFile, $outLog) = @_;
	open (my $in, "<", $fastaFile) or LOG($outLog, "Failed to open $fastaFile: $!\n") and exit 1;
	open (my $out, ">", "$fastaFile.out") or LOG($outLog, "Failed to write to $fastaFile.out: $!\n") and exit 1;
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$seq =~ tr/a-z/A-Z/;
		print $out "$def\n$seq\n";
	}
	close $in;
	close $out;
	system("/bin/mv $fastaFile.out $fastaFile") == 0 or LOG($outLog, "Failed to /bin/mv $fastaFile.out $fastaFile: $!\n") and exit 1;
	return("$fastaFile");
}

sub check_sanity {
	my ($opts) = @_;

	my ($usageshort, $usage, $usagelong, $example) = getUsage();

	print $example and exit if defined $opt_e;
	my $checkopts = 0;
	if (defined $opts) {
		my %opts = %{$opts};
		foreach my $key (sort keys %opts) {
			$checkopts ++ if defined $opts{$key};
		}
	}
	die "\n$usageshort\n" if $checkopts == 0;
	die "\n$usage\n"      if defined $opt_h;
	die "\n$usagelong\n"  if defined $opt_H;

	my $errors = "";
	my $read0 = defined($opt_r) ? $opt_r : "FALSE";

	# -r read not defined
	$errors .= "-r $LGN<read.fq>$N: is not defined.\n" if not defined($opt_r);
	$errors .= "-r $LGN<read.fq>$N: is defined ($opt_r) but file does not exist!\n" if defined $opt_r and not -e $opt_r;


	$errors .= "-n $LGN<output directory>$N is not defined.\n" if not defined($opt_n);
	my $message = "";
	if (defined $opt_n and not -d $opt_n) {
		$message = `mkdir $opt_n 2>&1`;
		chomp($message);
	}
	$errors .= "-n $LGN<output directory>$N is defined ($YW$opt_n$N) but folder does not exist and cannot be created.\n   -> mkdir $YW$opt_n$N returned $LCY$message$N.\n" if defined $opt_n and not -d $opt_n;
	$errors .= "-L $LGN<min read length (<integer>:in bp; <integer>p:% amplicon>$N must be positive integer!\n" if defined($opt_L) and ($opt_L =~ /^0+\.?0*[p]?$/ or $opt_L !~ /^\d+[p]?$/);

	$errors .= "-i $LGN<geneindex.bed>$N: is not defined.\n" if not defined($opt_i);
	$errors .= "-i <$LGN<geneindex.bed>$N: is defined ($YW$opt_i$N) but file does not exist!\n" if defined $opt_i and not -e ($opt_i);
	$errors .= "-g $LGN<ref_genome.fa>$N: is not defined.\n" if not defined($opt_g);
	$errors .= "-g $LGN<ref_genome.fa>$N: is defined ($YW$opt_g$N) but does not exist!\n" if defined $opt_g and not -e ($opt_g);

	# label error
	my $label;
	if (defined $opt_l) {
		$label = $opt_l;
	}
	elsif (defined $opt_n) {
		my ($readName) = getFilename($opt_n, "full");
		if ($readName =~ /PCB\d+_bcBC\d+_plasmid.+_+desc.+/i) {
			my ($PCB, $BC, $plasmid, $desc) = $readName =~ /PCB(\d+)_bcBC(\d+)_plasmid(.+)_desc([a-zA-Z0-9_]+)\.?/;
			if (defined $desc) {
				$label = "PCB$PCB\_bcBC$BC\_plasmid$plasmid\_desc$desc";
			}
		}
		elsif ($readName =~ /PCB\d+/i) {
			my ($PCB) = $opt_n =~ /PCB(\d+)/;
			if (defined $PCB) {
				$label = "PCB$PCB";
			}
		}
	}
	if (((defined $label and $label !~ /^PCB\d+(_bcBC\d+_plasmid.+_desc.+)$/) or not defined $label) and defined $opt_r) {
		my ($readName) = getFilename($opt_r, "full");
		if ($readName =~ /PCB\d+_bcBC\d+_plasmid.+_+desc.+/i) {
			my ($PCB, $BC, $plasmid, $desc) = $readName =~ /PCB(\d+)_bcBC(\d+)_plasmid(.+)_desc([a-zA-Z0-9_]+)\.?/;
			if (defined $desc) {
				$label = "PCB$PCB\_bcBC$BC\_plasmid$plasmid\_desc$desc";
			}
		}
		elsif ($readName =~ /PCB\d+/i) {
			my ($PCB) = $opt_r =~ /PCB(\d+)/;
			if (defined $PCB) {
				$label = "PCB$PCB";
			}
		}
	}
	$opt_l = $label if not defined $opt_l and defined $label;

	$errors .= "-l $LGN<label>$N: is not defined!\n\t-> $YW Has to be this exact format$N:.\n\t  - Invivo: use PCB<digits> e.g. PCB1 or PCB12636.\n\t  - Invitro: Use PCB<digits>_bcBC<digits>_plasmid<string>_desc<string> e.g. PCB13_bcBC13_plasmidPFC53_descSUPERCOILEDCIRC.\n" if not defined($opt_l);

	# label exists but wrong format
	$errors .= "-l $LGN<label$N: is defined ($YW$opt_l$N) but wrong format (${LPR}case sensitive$N)! PCB<digits> or PCB<digits>_bcBC<digits>_plasmid<string>_desc<string>.\n" if defined $opt_l and $opt_l !~ /^PCB\d+$/ and $opt_l !~ /^PCB\d+_bcBC\d+_plasmid[0-9a-zA-Z]+_desc[0-9a-zA-Z]+$/;

	if (not -d "$footLoopScriptsFolder/.sortTMP/") {
		system("mkdir $footLoopScriptsFolder/.sortTMP/") == 0 or $errors .= "Failed to make directory $footLoopScriptsFolder/.sortTMP/: $!\n";
	}
	
	if ($errors ne "") {
		print "$usageshort";
		print "\n${LRD}########## ${YW}FATAL ERROR${LRD} ##########\n\nREASON:$N\n";
		print $errors;
		print "\n\n${LRD}##################################$N\n\n" ;
		print "Use -Z if you're running invitro plasmid data!\n";
		die "\n";
	}
	print "Use -Z if you're running invitro plasmid data!\n";
}

sub record_options {
	my ($opts, $opts2, $outReadLog, $version, $outLog) = @_;
	my $optPrint = "$0";
	foreach my $opt (sort keys %{$opts}) {
		if (defined $opts->{$opt} and $opts->{$opt} eq "NONE") {
			$optPrint .= " -$opt" if defined $opts2->{$opt};
			LOG($outReadLog, "footLoop.pl,-$opt,\$TRUE\n","NA") if defined $opts2->{$opt};
			LOG($outReadLog, "footLoop.pl,-$opt,\$FALSE\n","NA") if not defined $opts2->{$opt};
		}
		elsif (defined $opts->{$opt} and $opts->{$opt} ne "NONE") {
			$optPrint .= " -$opt $opts->{$opt}";
			LOG($outReadLog, "footLoop.pl,-$opt,$opts->{$opt}\n","NA");
		}
	}
	my  $param = "
${YW}Initializing...$N

Version    : $version
Date       : $date
Run ID     : $uuid
Run script : $optPrint
	
${YW}Input Parameters$N
1. -r ${CY}Read/SAM$N  :$readFile
2. -g ${CY}Genome$N    :$genomeFile
3. -n ${CY}OutDir$N    :$outDir
4. -i ${CY}Index$N     :$geneIndexFile
5. -L ${CY}MinRdLen$N  :$minReadL
6. -q ${CY}MinMapQ$N   :$minMapQ
7. -Z ${CY}Bismark$N   :$bismarkOpt

";
	LOG($outLog, $param);
}

sub check_if_result_exist {
	my ($resFiles, $outLog) = @_;
	return if not defined $resFiles;
	my @resFiles = @{$resFiles};	
	return if @resFiles == 0;
	if (defined $opt_n and -d $opt_n and not defined $opt_F) {
		if (@{$resFiles} != 0) {
			my $resFileCheck = 0;
			foreach my $resFile (@{$resFiles}) {
				$resFileCheck = 1 if -e $resFile and -s $resFile > 0;
			}
			LOG($outLog, "-d $LGN$opt_n$N $resFiles->[0] already exist!\n") if $resFileCheck == 1;
			exit 0 if $resFileCheck == 1;
		}
		else {
			LOG($outLog, "End of run file does not exist!\n");
		}
	}
	return;
}


sub getUsage {
my $usageshort = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N [options..]
\t$CY-r$N read.fq
\t$LPR-n$N output_dir
\t$LGN-g$N genome.fa
\t$YW-i$N geneIndex.bed
\t$CY-x$N [0] left buffer in bp e.g. -x -10
\t$LPR-y$N [0] right buffer in bp e.g. -y 10
\t$LGN-L$N [95p] minimum read length (percent: 95p)
\t$YW-q$N [0] minimum map quality (use 0 since we're expecting a lot of conversions)
\t$CY-Z$N toggle to use non-stringent mapping
\t$LPR-F$N toggle to redo bismark mapping even if a .sam or .bam file is present in output_dir

Do $YW$0$N $LGN-h$N for longer explanation
Do $YW$0$N $LGN-e$N for example run

${LRD}IMPORTANT!!$N If you see a lot of 'Chromosomal sequence could not be extracted for..' try adding $YW-x -10 -y 10$N
-> If you still see then try adding $YW-x -50 -y 50$N

";

my $usage = "

$LGN-Z$N: Make mapping parameters not stringent. 
Bismark uses bowtie2 mapping parameters:
- normal stringency (default) formula: --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
- lowest stringency (-Z)      formula: --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8

-H: Longer explanation about bismark parameters above.

(..continued; Use -H to display longer help!)

";

my $usagelong = "$usage

${YW}Required:$N
1) -r: Input Fastq file
2) -n: Output directoy name (will create if not exist)
3) -g: Fasta file of reference genome (e.g. hg19.fa)
4) -i: index bed file (bed6 with gene name on each)

${YW}Optional [default]:$N

-x: ${LGN}[0]$N Add this number from the$GN start$N of the index (strand taken into account)

-y: ${LGN}[0]$N Add this number from the$GN end$N of the index, for example:
\tcoordinate is:$GN chr1 200 500 SEQ 0 +$N and$GN -x -100 -y 50$N becomes: chr1 100 550 SEQ 0 +$N
\tcoordinate is:$GN chr1 200 500 SEQ 0 -$N and$GN -x -100 -y 50$N becomes: chr1 150 600 SEQ 0 -$N

-p: (not implemented) Run bismark script in paralel.

-c: consider Cs in CpG context (default: off; Don't include Cs in CpG context)

-L: ${LGN}[95p]$N minimum ${CY}read$N length in bp to be considered valid read
    Add \"p\" to make it percent of amplicon length
    e.g. 95p = 95% of amplicon length.

-F: Force re-create BAM file

-q: ${LGN}[0]$N q = map quality, where probability of mapped incorrectly = 10^(q/-10). 
Some quick -q references:
- 0  = 100  %
- 10 = 10   %
- 20 = 1    %
- 30 = 0.1  %
- 40 = 0.01 %

${LCY}PS: Using -q 0 when mapping barcoded reads onto KNOWN unique amplicon sequences will not produce incorrect result$N.


FYI, Bowtie2 mapping parameters used in bismark2:

--rdg <int1>,<int2>      Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets a penalty
                         of <int1> + N * <int2>. Default: 5, 3.

--rfg <int1>,<int2>      Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of length N gets
                         a penalty of <int1> + N * <int2>. Default: 5, 3.

--score_min <func>       Sets a function governing the minimum alignment score needed for an alignment to be considered
                         \"valid\" (i.e. good enough to report). This is a function of read length. For instance, specifying
                         L,0,-0.2 sets the minimum-score function f to f(x) = 0 + -0.2 * x, where x is the read length.
                         See also: setting function options at http://bowtie-bio.sourceforge.net/bowtie2. The default is
                         L,0,-0.2.
(bowtie2)
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)

Stringent         : --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
Not stringent (-Z): --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8
--rdg 5,3: read gap of length N penalty = 5+3N (stringent) or 2+N (not stringent)
--rfg 5,3: ref. gap of length N penalty = 5+3N (stringent) or 2+N (not stringent)
--score_min: minimum threshold of a read with length of L: -0.3L (stringent) or -0.8L (not stringent)

";

my $example = "

${GN}Example:$N

## Download hg19.fa.gz from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz, gunzip it, then put into footLoop folder.
## Download example.tar.gz and move everything in it into footLoop folder.
tar zxvf example.tar.gz
mv example/* footLoop/
cd footLoop
./footLoop.pl -r PCB190425.fq.gz -n PCB190425_MAP -g hg19.fa -l PCB190425 -x -10 -y 10 -i geneIndexes.bed
./footPeak.pl -n PCB190425_MAP -o PCB190425_PEAK
./footClust.pl -n PCB190425_PEAK
./footPeak_graph.pl -n PCB190425_PEAK
./footPeakGTF.pl -n PCB190425_PEAK
./footPeak_GCprofile.pl -n PCB190425_PEAK -i geneIndexes.bed
./footStats.pl -n PCB190425_PEAK

";

	return($usageshort, $usage, $usagelong, $example);
}
