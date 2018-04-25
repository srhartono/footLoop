#!/usr/bin/perl
# Version 160831_Fixed_PrintOutput at the same file (step 8)
use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
#use Getopt::Std::WithCheck;
use vars   qw($opt_v $opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f $opt_l);
my @opts = qw($opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_l);
getopts("vr:g:i:n:L:x:y:q:HhZFpl:");

BEGIN {
	my ($bedtools) = `bedtools --version`;
	my ($bowtie2) = `bowtie2 --version`;
	my ($bismark) = `bismark --version`;
	my ($bismark_genome_preparation) = `bismark_genome_preparation --version`;
	
	if (not defined $bedtools or $bedtools =~ /command not found/ or $bedtools =~ /bedtools v?([01].\d+|2\.0[0-9]|2\.1[0-6])/) {
		print "Please install bedtools at least version 2.17 before proceeding!\n";
		$bedtools = 0;
	}
	if (not defined $bowtie2 or $bowtie2 =~ /command not found/ or $bowtie2 =~ /version [0-1]./) {
		print "Please install bowtie2 at least version 2.1.0 before proceeding!\n";
		$bowtie2 = 0;
	}
	if (not defined $bismark or $bismark =~ /command not found/ or $bismark =~ /v?(0\.1[0-2]|0\.0[0-9])/) {
		print "Please install bismark at least version 0.13 before proceeding!\n";
		$bismark = 0;
	}
	if (not defined $bismark_genome_preparation or $bismark_genome_preparation =~ /command not found/ or $bismark_genome_preparation =~ /v?(0\.1[0-2]|0\.0[0-9])/) {
		print "\n\nPlease install bismark_genome_preparation at least version 0.13 before proceeding!\n\n";
		$bismark_genome_preparation = 0;
	}
	print "- bedtools v2.17+ exists:" . `which bedtools` if $bedtools ne 0;
	print "- bowtie2 v2.1+ exists:" . `which bowtie2` if $bowtie2 ne 0;
	print "- bismark v0.13+ exists:" . `which bismark` if $bismark ne 0;
	print "- bismark_genome_preparation v0.13+ exists:" . `which bismark_genome_preparation` if $bismark_genome_preparation ne 0;
	die if $bedtools eq 0 or $bowtie2 eq 0 or $bismark eq 0 or $bismark_genome_preparation eq 0;
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";
($opt_r, $opt_i, $opt_n, $opt_g, $opt_x, $opt_y) = run_example() if @ARGV and $ARGV[0] eq "ex";

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

###################
# 0. Check Sanity #
##################

my %OPTS  = ('r' => $opt_r, 'g' => $opt_g, 'i' => $opt_i, 'n' => $opt_n, 
				 'L' => $opt_L, 'x' => $opt_x, 'y' => $opt_y, 'p' => $opt_p, 
				 'q' => $opt_q, 'F' => 'NONE', 'Z' => 'NONE',
				 'p' => 'NONE', 'q' => $opt_q, 'l' => $opt_l);
my %OPTS2 = ('p' => $opt_p, 'Z' => $opt_Z, 'F' => $opt_F);

sanityCheck(\%OPTS, \%OPTS2);

###################
# 1. Define Input #
###################
my ($readFile, $genomeFile, $geneIndexFile, $outDir) = getFullpathAll($opt_r, $opt_g, $opt_i, $opt_n);
my ($readFilename)  = getFilename($opt_r, "full");
my ($geneIndexName) = getFilename($geneIndexFile);
my $minMapQ    = (not defined($opt_q)) ? 0  : $opt_q;
$opt_L = "85p" if not defined $opt_L;
my ($minReadL)   = (not defined($opt_L)) ? 85 : $opt_L =~ /p$/i ? $opt_L =~ /^(.+)p$/i : $opt_L;
my $bufferL    = (not defined($opt_x)) ? 0  : $opt_x;
my $bufferR    = (not defined($opt_y)) ? 0  : $opt_y;
my $readName   = getFilename($readFile, "full");
my $bismarkOpt = (defined $opt_Z) ? "--bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8" : "--bowtie2 --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3";
my $samFile	   = ($readFile =~ /.f(ast)?q(.gz)?$/) ? $outDir .  "/$readName\_bismark_bt2.sam" : $readFile;
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
   die "Please make sure your output folder (-o) contain PCB(number) e.g. PCB12: 180202_PCB12_footloop_output\n\n";
}

if (not defined $opt_l and $outDir =~ /_BC\d+_PFC\d+_\w+\./i) {
   my ($label2) = $outDir =~ /_(BC\w+)\./i;
   $label = "$label\_$label2";
   $label = uc($label);
   $label =~ s/PCB[0\-_]+(\d+)/PCB$1/g;
}
	$OPTS{l} = $label;

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

# Add buffers to geneIndexFile (bedtools_bed_change.pl) and get their sequences using fastaFromBed, then get their info from fasta=$seqFile

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

# Do footLoop_2_sam_to_peak.pl
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
	               png(PNGout[i],height=4*2*dim(df)[1],width=2*dim(df)[2])
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
	#	}
		close $out;
		my $cmd = "R --no-save < $Rscript >> $logFile\n";
		print ($cmd);
		system($cmd) if -e $Rscript and -s $Rscript > 10;
		system("rm Rplots.pdf") if -e "Rplots.pdf";
		LOG($outLog, "\n${LPR}If there is not PDF made from step (8) then it's due to too low number of read/peak$N\n\n");
	}
}

sub run_example {
	my $exFolder = "$footLoopDir/example";
	print "Running example folder $LCY$exFolder$N\n";
	my $opt_r = "$exFolder/CALM3/CALM3.sam";
	my $opt_i = "$exFolder/CALM3/CALM3index.bed";
	my $opt_g = "$exFolder/CALM3/hg19.fa";
	my $opt_n = "$exFolder/CALM3/CALM3out/";
	return($opt_r, $opt_i, $opt_n, $opt_g, -100, 100);
}

sub fix_samFile {
	my ($samFile, $seqFile, $outReadLog, $outLog) = @_;
	LOG($outLog, "\n\ta. Fixing sam file $CY$samFile$N with footLoop_2_sam_to_peak.pl\n");
	my ($samMD5) = getMD5($samFile);
	my $origDir = "$outDir/.0_orig_$samMD5/";
	check_if_result_exist(["$origDir/.GOOD"], $outLog);
	makedir("$origDir") if not -d "$origDir";
	my $checkSam;
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
			$checkSam = $samLineCount1 - 500 > $samLineCount2 ? 0 : 2;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: fixed sam file $LCY$origDir/$samFileName.fixed.gz$N exists but total row is less than total samFile $samFile row ($samLineCount1 - 500 > samFile.fixed.gz: $samLineCount2)!\n") if $checkSam == 0;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile::$LGN SUCCESS!!$N fixed sam file $LCY$origDir/$samFileName.fixed.gz$N exists (MD5=$LGN$samMD5$N) and total row $LGN($samLineCount2)$N >= total samFile row $LGN($samLineCount1 - 500)$N ($LCY$samFile$N)!\n") if $checkSam == 2;
		}
		if (-e "$origDir/$samFileName.fixed" and $checkSam == 0) {
			my ($samLineCount2) = linecount("$origDir/$samFileName.fixed");
			my ($samLineCount1) = linecount($samFile);
			$checkSam = $samLineCount1 - 500 > $samLineCount2 ? 0 : 1;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: .gz does not exist and fixed sam file $LCY$origDir/$samFileName.fixed$N exists but total row is less than total samFile $samFile row ($samLineCount1 - 500 > samFile.fixed: $samLineCount2)!\n") if $checkSam == 0;
			LOG($outLog, "\tfootLoop.pl subroutine fix_samFile::$LGN SUCCESS!!$N fixed sam file $LCY$origDir/$samFileName.fixed$N exists (MD5=$LGN$samMD5$N) and total row $LGN($samLineCount2)$N >= total samFile row $LGN($samLineCount1 - 500)$N ($LCY$samFile$N)!\n") if $checkSam == 1;
		}
	}
	if ($checkSam == 0) {
		LOG($outLog, "\tfootLoop.pl subroutine fix_samFile:: fixed sam file $LCY$origDir/$samFileName.fixed$N or .gz does not exist!\n");
		LOG($outLog, "\t${YW}footLoop_2_sam_to_peak.pl -f $outDir -s $seqFile -o $origDir$N\n");
		system("footLoop_2_sam_to_peak.pl -f $outDir -o $origDir") == 0 or LOG($outLog, "Failed to run footLoop_2_sam_to_peak.pl -f $outDir -o $origDir: $!\n") and exit 1;
		LOG($outReadLog, "footLoop.pl,fix_samFile,footLoop_2_sam_to_peak.pl -f $outDir -o $origDir\n","NA");
		LOG($outLog, "\tgzip $origDir/$samFileName.fixed");
		system("gzip -f $origDir/$samFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $origDir/$samFileName.fixed: $!\n");
		$checkSam = 1;
	}
	else {
		LOG($outLog, "\t${LGN}WAS NOT RUN$N: ${YW}::: footLoop_2_sam_to_peak.pl -f $outDir -s $seqFile -o $origDir :::$N\n");
		# rm old (bad) .gz if it exists
		LOG($outLog, "\t/bin/rm $samFileGZ") if -e '$samFileGZ';
		system("/bin/rm $samFileGZ") == 0 or LOG($outLog, 'Failed to rm $samFileGZ: $!\n') if -e '$samFileGZ';
		# gzip the new .fixed
		LOG($outLog, "\tgzip $origDir/$samFileName.fixed");
		system("gzip -f $origDir/$samFileName.fixed") == 0 or LOG($outLog, "\tFailed to gzip $origDir/$samFileName.fixed: $!\n");
	}

	# re-md5 samfile gz
	LOG($outLog, "\t${YW}md5sum $samFileGZ > $origDir/.$samFileName.fixed.gz.md5$N\n");
	system("md5sum $samFileGZ > $origDir/.$samFileName.fixed.gz.md5") == 0 or LOG($outLog, "Failed to md5sum $samFileGZ > $origDir/.$samFileName.fixed.gz.md5: $!\n") and exit 1;
	LOG($outReadLog, "footLoop.pl,samFixedFile,$samFileGZ\n","NA");
	LOG($outReadLog, "footLoop.pl,samFixedFileMD5,$samMD5\n","NA");
	return($samFile, $origDir);
}

sub parse_samFile {
	my ($samFile, $seqFile, $outReadLog, $outLog) = @_;
	my ($samFileName) = getFilename($samFile);

	LOG($outLog, "\ta. Parsing sam file $CY$samFile$N and getting only high quality reads\n");
	open(my $notused, ">", "$outDir/.$readFilename.notused") or LOG($outLog, "Cannot open $outDir/.$readFilename.notused: $!\n") and exit 1;
	open(my $sam, $samFile) or LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!") and exit 1;
	
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
			print $notused "\t$CY$read$N length of seq ($CY" . $seqLen. "$N) is less than $SEQ->{$gene}{minReadL} bp (length of original sequence is ($CY" . $SEQ->{$gene}{geneL} . "$N)!\n";
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
	open (my $inSamFix, "zcat $origDir/$samFileName.fixed.gz|") or LOG($outLog, "Failed to open $origDir/$samFileName.fixed.gz: $!\n") and exit 1;
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
		if ($type eq "6_BOTH") {
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
		#print "NUC=$nuc\n";die;
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
#			print "\t- $nuc\t$beg\t$end\t$prev0$YW$curr$N$next0\n";
		}
	}
	$data =~ s/^,// if defined $data;
	return($data);
}

sub make_bismark_index {
	my ($geneIndexFa, $bismark_folder, $bismarkOpt, $outLog) = @_;
	LOG($outLog, "\n\ta. Running bismark_genome_preparation$N --bowtie2 $bismark_folder$N\n");
	my $run_boolean = "\t${LGN}WAS NOT RUN$N:${YW} ";
	my $cmd = "bismark_genome_preparation --bowtie2 $bismark_folder > $bismark_folder/LOG.txt 2>&1 && md5sum $geneIndexFa > $bismark_folder/Bisulfite_Genome/.MD5SUM";
#"bismark_genome_preparation --bowtie2 $bismark_folder && md5sum $geneIndexFa > $bismark_folder/Bisulfite_Genome/.MD5SUM";
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
#		die "bismark folder = $LCY$bismark_folder$N\n";
		LOG($outLog, "\tEither bismark folder didn't exist or older bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n") if defined $md5sum;
		system($cmd) == 0 or die "Failed to run bismark genome preparation: $!\n";
	}
	else { 
		LOG($outLog, "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n");
	}
	LOG($outLog, "${run_boolean} $YW ::: $cmd :::$N\n");
#	die "OK\n";#bismark folder = $LCY$bismark_folder$N\n";
}

sub run_bismark {
	my ($readFile, $outDir, $mysam, $opt_F, $outReadLog, $outLog) = @_;
	my $MAP = "";
	LOG($outLog, "\n\tb. Running bismark\n");
	my ($mysamFilename) = getFilename($mysam, "full");
	my $run_boolean = "\n\t${LGN}WAS NOT RUN$N:${YW} ";
	if (-e $mysam and not -e "$outDir/$mysamFilename") {
		system("/bin/ln -s $mysam $outDir/$mysamFilename") == 0 or LOG($outLog, "Failed to /bin/ln $mysam $outDir/$mysamFilename: $!\n") and exit 1;
	}
	if (defined $opt_F or not -e $mysam) {
		$run_boolean = "\n$YW\t ";
		if ($opt_p) {
			my $outFolder = $outDir . "/.bismark_paralel/";
			if (-d $outFolder) {
				my @files = <$outFolder/*>;
				if (@files != 0) {
#					my $date2 = date();
#					makedir("$outFolder/backup_$date2/");
#					LOG($outLog, "\tmoving all in $outFolder into $outFolder/backup (with date=$date2)\n\t$YW/bin/mv $outFolder/* $outFolder/backup_$date2/$N");
					system("/bin/rm $outFolder/*");# $outFolder/backup_$date2/");
				}
			}
			makedir($outFolder) if not -d $outFolder;
			LOG($outLog, "\t  Splitting $CY$readFile$N by 1000 sequences!\n\n####### SPLITRESULT LOG ########");	
			my $splitresult = `SplitFastq.pl -i $readFile -o $outFolder -n 1000`;
			LOG($outReadLog, "footLoop.pl,run_bismark,SplitFastq.pl -i $readFile -o $outFolder -n 1000\n");

			LOG($outLog, "$splitresult\n");	
			print "####### SPLITRESULT LOG #######\n\n\t  Running bismark in paralel!\n";
			my $result = system("run_script_in_paralel2.pl -v \"srun -p high --mem 8000 bismark -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\" $outFolder .part 20");
			LOG($outReadLog, "footLoop.pl,run_bismark,\"run_script_in_paralel2.pl -v \\\"srun -p high --mem 8000 bismark -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\\\" $outFolder .part 20");
			my @partSam = <$outFolder/*.part_bismark*.sam>; my $totalPartSam = @partSam;
			LOG($outLog, "\t  All part.sam has been made (total = $totalPartSam). Now making $CY$mysam$N and then removing the part sam\n");
			my @HEADER; my @REPORT;
			for (my $p = 0; $p < @partSam; $p++) {
				my $partSam = $partSam[$p];
				print "\t\tPutting $partSam into $mysam and removing it!\n";
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p == 0;
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p != 0;
				LOG($outReadLog, "footLoop.pl,run_bismark,cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mysam") if $p == 0;
				LOG($outReadLog, "footLoop.pl,run_bismark,cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mysam") if $p != 0;
				my ($bismark_report) = $partSam =~ /^(.+).sam/; $bismark_report .= "_SE_report.txt";
				my ($header, $report) = parse_bismark_report($bismark_report);
				my @header = @{$header};
				my @report = @{$report};
				@HEADER = @header if $p == 0;
				for (my $q = 0; $q < @header; $q++) {
					die "undefined $q header\n" if not defined $header[$q];
					die if $header[$q] ne $HEADER[$q];
					$REPORT[$q] += $report[$q] if $header[$q] !~ /^perc_/;
					$REPORT[$q] += $report[0]*$report[$q] if $header[$q] =~ /^perc_/;
				}
#				print "\t- Removing $CY$partSam$N: /bin/rm $partSam\n";
#				system("/bin/rm $partSam") == 0 or die "Failed to /bin/rm $partSam: $!\n";
			}
			for (my $q = 0; $q < @HEADER; $q++) {
				$REPORT[$q] = int(100*$REPORT[$q]/$REPORT[0]+0.5)/100 if $HEADER[$q] =~ /^perc_/;
			}
			$MAP  = "footLoop.pl,map," . "header\tlabel\t" . join("\t", @HEADER) . "\tfootLoop_outDir\tuuid\n";
		   $MAP .= "footLoop.pl,map," . "record\t$label\t" . join("\t", @REPORT) . "\t$outDir\t$uuid\n";
		}
		else {
			LOG($outLog, "\t  bismark -o <outDir> $LCY$bismarkOpt$N <bismark_folder> <readFile> > <outDir/.bismark_log> 2>&1\n");
			LOG($outReadLog, "footLoop.pl,bismark,bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1\n","NA");
			my $result = system("bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");

			if ($result != 0) {
				LOG($outLog, "\t\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_folder\n");
				make_bismark_index($seqFile, $bismark_folder, $bismarkOpt, $outLog);
				system("bismark_genome_preparation --bowtie2 $bismark_folder") == 0 or die "Failed to run bismark genome preparation: $!\n";
				LOG($outReadLog, "footLoop.pl,bismark,bismark_genome_preparation --bowtie2 $bismark_folder","NA");
				system("bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
				LOG($outReadLog, "footLoop.pl,bismark,bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");
			}
			LOG($outLog, "\t${GN}SUCCESS$N: Output $mysam\n");
			my ($bismark_report) = "$mysam" =~ /^(.+).sam/; $bismark_report .= "_SE_report.txt";
			my ($header, $report, $MAPTEMP) = parse_bismark_report($bismark_report);
			$MAP = $MAPTEMP;
		}
	}
	else {
		print "A\n";
		LOG($outLog, "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n");
		my ($bismark_report) = "$mysam" =~ /^(.+).sam/; $bismark_report .= "_SE_report.txt";
		my ($header, $report, $MAPTEMP) = parse_bismark_report($bismark_report);
		$MAP = $MAPTEMP;
	}
	LOG($outLog, "${run_boolean}::: bismark $bismarkOpt $bismark_folder $readFile :::$N\n");

	# print to $outReadLog
	LOG($outReadLog, "footLoop.pl,samFile,$outDir/$mysamFilename\n","NA");
	LOG($outReadLog, $MAP,"NA");

	return("$outDir/$mysamFilename");
}

sub parse_bismark_report {
	my ($bismark_report_file) = @_;
	open (my $inz, "<", $bismark_report_file) or LOG($outLog, "Failed to read from $LCY$bismark_report_file$N: $!\n") and die;
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
#	print join("\t", @header) . "\n" . join("\t", @report) . "\n\n";
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
		LOG($outLog, "$YW\t::: bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1 :::$N\n") == 0 or LOG($outLog, "Failed to get (beg $bufferL, end $bufferR) bp of $geneIndexFile!\n") and exit 1;
		system("bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1") == 0 or LOG($outLog, "\tfootLoop.pl::get_geneIndex_fasta: Failed to$YW bedtools_bed_change.pl -m -x $bufferL -y $bufferR -i $geneIndexFile -o $geneIndexFileNew >> $logFile 2>&1$N\n: $!\n") and exit 1;
		LOG($outLog, "\t${GN}SUCCESS$N: Created index file $LGN$geneIndexFileNew$N from bedtools bed change of $LCY$geneIndexFile$N\n");
	}
	return($geneIndexFileNew);
}

sub parse_geneIndexFile {
	my ($geneIndexFile, $genomeFile, $outDir, $minReadL, $outReadLog, $outLog) = @_;

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

#=SUB convert_seq {

sub sanityCheck {
	print "\n\nUse -Z if you're running invitro plasmid data!\n\n";
	my ($opts) = @_;

	my ($usageshort, $usage, $usagelong) = getUsage();
	my $checkopts = 0;
	foreach my $key (sort keys %{$opts}) {
		$checkopts ++ if defined $opts->{$key};
	}
	die "\n$usageshort\n" if $checkopts == 0;
	die "\n$usage\n"      if defined $opt_h;
	die "\n$usagelong\n"  if defined $opt_H;

	my $read0 = defined($opt_r) ? $opt_r : "FALSE";

	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -r <read.fq> not defined\n\n############################################\n\n" if not defined($opt_r);
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -r $opt_r DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if defined $opt_r and not -e $opt_r;
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -n <output directory> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_n);
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -i <geneindex.bed> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_i);
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -L <min read length in bp> must be positive integer!\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if defined($opt_L) and $opt_L !~ /p$/ and $opt_L !~ /^\d+$/;
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -g <ref_genome.fa [hg19.fa]> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_g);
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -i $opt_i DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not -e ($opt_i);
	die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -g $opt_g DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not -e ($opt_g);
	if (not -d "$footLoopDir/.sortTMP/") {
		system("mkdir $footLoopDir/.sortTMP/") == 0 or die "Failed to make directory $footLoopDir/.sortTMP/: $!\n";
	}
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
my $usageshort = "\n$LGN----------------------------$N\n
Usage: $YW$0$N [options..]
\t$CY-r$N read.fq
\t$LPR-n$N output_dir
\t$LGN-g$N genome.fa
\t$YW-i$N geneIndex.bed
\t$CY-x$N [0] left buffer in bp e.g. -100 for CALM3
\t$LPR-y$N [0] right buffer in bp e.g. 100 for CALM3
\t$LGN-L$N [500] minimum read length (percent: 50p)
\t$YW-q$N [0] minimum map quality (0-50)
\t$CY-Z$N toggle to use non-stringent mapping (-H for more explanation)
\t$LPR-F$N toggle to redo bismark mapping even if a .sam file is present in output_dir

Do $YW$0$N $LGN-h$N for longer explanation

${LRD}IMPORTANT!!$N If you see a lot of 'Chromosomal sequence could not be extracted for..' try adding $YW-x -10 -y 10$N
If you still see then try adding $YW-x -50 -y 50$N

$LGN----------------------------$N\n
";

#Usage: $YW$0$N ${GN}[options]$N [-r $CY<read.fq>$N OR -S $CY<read.bismark.sam>$N -n $PR<output folder>$N -g $GN<genome.fa>$N -i $YW<geneCoordinates.bed>$N
my $usage = "
${LGN}Options [default]:$N
-x $LGN<Left pad [0]>$N -y $LGN<Right pad [0]>$N -L $LGN<min READ length [500]>$N -q $LGN<min map quality [0]>$N

Example:
$YW$0$N\n\t$CY-r$N example/pacbio12ccs.fq \n\t$PR-n$N myoutput\n\t$GN-g$N example/hg19.fa.fa\n\t$YW-i$N example/pacbio12index.bed\n\t-x-100\n\t-y 100\n\t-L 1000\n

$LGN-Z$N: Make mapping parameters not stringent. 
Bismark uses bowtie2 mapping parameters:
- normal stringency (default) formula: --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
- lowest stringency (-Z)      formula: --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8

-H: Longer explanation about bismark parameters above.


(..continued; Use -H to display longer help!)

";

my $usagelong = "$usage

${YW}Required:$N
1) Either:
	${CY}-r$N: must be path to reads file itself named pacbio.fastq
	${CY}-S$N: supply a sam file
2) -n: Output directoy name (will create if not exist)
3) -g: Fasta file of reference genome (e.g. hg19.fa)
4) -i: index bed file (bed4 with gene name on each)

${YW}Optional [default]:$N
-x: ${LGN}[0]$N Add this number from the$GN start$N of the index (strand taken into account)
-y: ${LGN}[0]$N Add this number from the$GN end$N of the index 
e.g. seq is:$GN chr1 200 500 SEQ 0 +$N and$GN -x -100 -y 50$N becomes: chr1 100 550 SEQ 0 +$N
e.g. seq is:$GN chr1 200 500 SEQ 0 -$N and$GN -x -100 -y 50$N becomes: chr1 150 600 SEQ 0 -$N
-p: Run bismark script in paralel
-H: Run HMM peak caller$LRD NOT IMPLEMENTED$N
-c: <default: don't include CpG> consider Cs in CpG context
-t: ${LGN}[0.65]$N percentage (0.0-1.0) of Cs that must be converted to be considered an R-loop
-l: ${LGN}[100]$N minimum length in base pairs to be considered an R-loop peak (also for HMM)
-L: ${LGN}[500]$N minimum ${CY}read$N length in bp to be considered valid read
    Add \"p\" to make it percent of amplicon length
    e.g. 50p = 50% of amplicon length.

If you want to re-use previously made files from step 4 and/or 6, use:
-1: For both step 4 and step 6 (Positive/NegativeFinal.txt and methylationPos/Neg.txt)
-4: For step 4 (Positive/NegativeFinal.txt)
-6: For step 6 (methylationPos/Neg.txt)
$LRD IF -1 is present, it'll override -4 and -6!
-F: Force re-create SAM file
-q: ${LGN}[0]$N q = map quality, where probability of (wrong) = 10^(q/-10). 
However for something like R-loop footprint from pacbio, where most stuff will be weird (indel+converted), this score is kind of meaningless so use default (0)
Some examples:
-  0 = 100     %
- 10 =  10     %
- 20 =   1     %
- 30 =    0.1  %
- 40 =    0.01 %

${GN}Example:$N

If you have .fq file but no SAM file:
$CY$0$N ${YW}-r example/pacbio12ccs.fq$N -n ${CY}myoutput$N -g$CY example/hg19.fa.fa$N -i$CY example/pacbio12index.bed$N -p -t ${CY}0.65$N -L$CY 500$N -l$CY 250$N 

If you have SAM file, use -S instead of -r (everything is the same as above except$YW yellow$N):
$CY$0$N ${YW}-S example/pacbio12ccs.bismark.sam$N -n ${CY}myoutput$N -g$CY example/hg19.fa.fa$N -i$CY example/pacbio12index.bed$N -p -t ${CY}0.65$N -L$CY 500$N -l$CY 250$N 



###### BISMARK (BOWTIE2) PARAMETERS
(bowtie2 in bismark)
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

pFC53_small.bed (from PCB14) is 2500bp long, with 531 C. 
All converted Cs and a 5\% mismatches (0.05*(2500-531)) and 1\% gap= 100bp mismatches
";

return($usageshort, $usage, $usagelong);
}
__END__

	if (defined $opt_1 or defined $opt_4) {
		my $dir4 = defined $opt_1 ? $opt_1 : $opt_4;
		my $PositiveFinal  = $dir4 . "/PositiveFinal.txt";
		my $NegativeFinal  = $dir4 . "/NegativeFinal.txt";
		die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N" . "-1 or -4 is given ($dir4) but $PositiveFinal does not exist!\n" if not -e $PositiveFinal;
		die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N" . "-1 or -4 is given ($dir4) but $NegativeFinal does not exist!\n" if not -e $NegativeFinal;
	}
	if (defined $opt_1 or defined $opt_6) {
		my $dir6 = defined $opt_1 ? $opt_1 : $opt_6;
		my $MethylationPos = $dir6 . "/methylationPos.txt";
		my $MethylationNeg = $dir6 . "/methylationNeg.txt";
		die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N" . "-1 or -6 is given ($dir6) but $MethylationPos does not exist!\n" if not -e $MethylationPos;
		die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N" . "-1 or -6 is given ($dir6) but $MethylationNeg does not exist!\n" if not -e $MethylationNeg;		
	}
}

__END__
Nov2017

					pdf(paste(files[i],\"_$threshold.pdf\",sep=\"\"))
					heatmap.3(
					x=df,
					dendrogram=\"none\",
					Rowv=myclust, Colv=FALSE,
					labRow=TRUE,labCol=FALSE,
					sideRow=2,
					cexRow=0.2,
					ColIndividualColors=seq,
					breaks=$breaks,
					color.FUN=function(x) $colors,
					main=files[i],
					sub=paste(\"(\",dim(df)[1],\" peaks)\",sep=\"\"),
					cex.main = 0.5
					)
					dev.off()

					print(paste(\"   \",i,\"B. Doing PNG of \",files[i],sep=\"\"))
					png(paste(files[i],\"_$threshold.png\",sep=\"\"),height=2000,width=2000)
					heatmap.3(
					x=df,
					dendrogram=\"none\",
					Rowv=myclust, Colv=FALSE,
					labRow=TRUE,labCol=FALSE,
					sideRow=2,
					cexRow=0.2,
					ColIndividualColors=seq,
					breaks=$breaks,
					color.FUN=function(x) $colors,
					main=files[i],
					sub=paste(\"(\",dim(df)[1],\" peaks)\",sep=\"\"),
					cex.main = 0.5
					)

					dev.off()


__END__
#CUT1
#AIRN_PFC53_REVERSE      1275-1287
#AIRN_PFC53_REVERSE      1372-1382
#AIRN_PFC53_REVERSE      1490-1507

						#print "CHUNK\t$chunk\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
						my ($begPeak, $begPeak2) = $chunk =~ /^([023678]*)(1|4|5|9)/;
						#print "Begpeak\t$begPeak\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
						$begPeak = defined $begPeak ? length($begPeak) + $i + 1 : $i + 1;
						my ($temp1, $temp2, $midPeak, $temp3, $temp4) = $chunk =~ /^([023678]*)(1|4|5|9)(.+)(1|4|5|9)([023678]*)$/;
						$midPeak = "undefined!\n" if not defined $midPeak;
						#print "Midpeak\t$midPeak\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
						my $printPeak = "$midPeak";
						my ($endPeak, $endPeak2, $endPeak3) = $chunk =~ /^(.*)(1|4|5|9)([023678]*)$/;
						#print "endpeak\t$endPeak\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
						$endPeak = defined $endPeak ? length($endPeak) + $i + 1 : $i + 1;
						my $checkPeak = 0;
						#print "\tGENE=$gene has peak! (conv = $conPer > threshold=$opt_t), checking at beg=$LGN$begPeak$N, end=$LCY$endPeak$N\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
						foreach my $beg_C (sort keys %{$lotsOfC{$gene}}) {
							my $end_C = $lotsOfC{$gene}{$beg_C};
							if (($begPeak >= $beg_C and $begPeak <= $end_C) or ($endPeak >= $beg_C and $endPeak <= $end_C)) {
								$printPeak .= "\tgene=$gene begPeak=$begPeak endPeak=$endPeak is within $beg_C-$end_C\n";# if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
								$checkPeak = 1; last;
							}
						}
#						die if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
#CUT1

__END__
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -n <output directory> not defined\n\n$usageshort\n\n" if not defined($opt_n);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r <read.fq> and $opt_s <read.sam> both not defined\n\n" if not defined($opt_r) and not defined($opt_S);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i <geneindex.bed> not defined\n\n" if not defined($opt_i);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -l <min peak length in bp> must be positive integer!\n\n" if defined($opt_l) and $opt_l !~ /^\d+$/;
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -L <min read length in bp> must be positive integer!\n\n" if defined($opt_L) and $opt_L !~ /p$/ and $opt_L !~ /^\d+$/;
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g <ref_genome.fa [hg19.fa]> not defined\n\n" if not defined($opt_g);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r $read0 and -S $sam0 both DOES NOT EXIST (provide at least one!)\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if (defined($opt_r) and not -e ($opt_r)) or (defined($opt_S) and not -e ($opt_S));
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i $opt_i DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_i);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g $opt_g DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_g);
#die "\n########## USAGE ##########\n$N$usage \n${LRD}########## FATAL ERROR ##########\n\n!!!$N -S $CY$opt_S$N exists but seems to be empty!\n$N$usage" if defined($opt_S) and (not -e $opt_S or (-s $opt_S < 10));



__END__
HERE1
	if (defined $opt_1 or defined $opt_4) {
		my $finalPositive = $outDir . "/PositiveFinal.txt";
		my $finalNegative = $outDir . "/NegativeFinal.txt";
		print "\tWarning: Step 4: either -1 or -4 is given, but $finalPositive and $finalNegative already exists! Moved these into $finalPositive.backup and $finalNegative.backup!\n" if (-e $finalPositive or -e $finalNegative);
		system("/bin/mv $finalPositive $finalPositive.backup") == 0 or die "Failed to mv $finalPositive to $finalPositive.backup: $!\n" if -e "$finalPositive";
		system("/bin/mv $finalNegative $finalNegative.backup") == 0 or die "Failed to mv $finalNegative to $finalNegative.backup: $!\n" if -e "$finalNegative";	

		# Create symbolic link
		my $dir4 = defined $opt_1 ? $opt_1 : $opt_4;
		my $finalPositive2  = $dir4 . "/PositiveFinal.txt";
		my $finalNegative2  = $dir4 . "/NegativeFinal.txt";
		($finalPositive2) = myFootLib::getFullpath($finalPositive2);
		($finalNegative2) = myFootLib::getFullpath($finalNegative2);
	
		system("ln -s $finalPositive2 $finalPositive") == 0 or die "\tFailed at step 4: Cannot create link of $finalPositive2 to $outDir: $!\n";
		system("ln -s $finalNegative2 $finalNegative") == 0 or die "\tFailed at step 4: Cannot create link of $finalNegative2 to $outDir: $!\n";
	}

	if (defined $opt_1 or defined $opt_6) {
		# If in new folder, finalPositive already exist, then die UNLESS -f is given
		my $methylationPos = defined $opt_c ? $outDir . "/methylatioPosCG.txt" : $outDir . "/methylationPos.txt";
		my $methylationNeg = defined $opt_c ? $outDir . "/methylatioNegCG.txt" : $outDir . "/methylationNeg.txt";
		print "\tWarning: Step 6: either -1 or -6 is given, but $methylationPos and $methylationNeg already exists! Moved these into $methylationPos.backup and $methylationNeg.backup!\n" if (-e $methylationPos or -e $methylationNeg);
		system("/bin/mv $methylationPos $methylationPos.backup") == 0 or die "Failed to move $methylationPos to $methylationPos.backup: $!\n" if -e "$methylationPos";
		system("/bin/mv $methylationNeg $methylationNeg.backup") == 0 or die "Failed to move $methylationNeg to $methylationNeg.backup: $!\n" if -e "$methylationNeg";	

		# Create symbolic link
		my $dir6 = defined $opt_1 ? $opt_1 : $opt_6;
		my $methylationPos2 = defined ($opt_c) ? $dir6 . "/methylationPosCG.txt" : $dir6 . "/methylationPos.txt";
		my $methylationNeg2 = defined ($opt_c) ? $dir6 . "/methylationNegCG.txt" : $dir6 . "/methylationNeg.txt";
		($methylationPos2) = myFootLib::getFullpath($methylationPos2);
		($methylationNeg2) = myFootLib::getFullpath($methylationNeg2);
		system("ln -s $methylationPos2 $outDir/") == 0 or die "\tFailed at step 6: Cannot create link of $methylationPos2 to $outDir: $!\n";
		system("ln -s $methylationNeg2 $outDir/") == 0 or die "\tFailed at step 6: Cannot create link of $methylationNeg2 to $outDir: $!\n";
	}
HERE1END

PART6
	#my $finalPos = $outDir . "$gene\_Pos" . ($opt_t*100) . ".txt";
	#my $finalNeg = $outDir . "$gene\_Neg" . ($opt_t*100) . ".txt";
	#$finalPos = $outDir . "$gene\_Pos" . ($opt_t*100) . "CG.txt" if($opt_c);
	#$finalNeg = $outDir . "$gene\_Neg" . ($opt_t*100) . "CG.txt" if($opt_c);
	#my $finalPos_NOPEAK = $outDir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK.txt";
	#my $finalNeg_NOPEAK = $outDir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK.txt";
	#$finalPos_NOPEAK = $outDir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK_CG.txt" if($opt_c);
	#$finalNeg_NOPEAK = $outDir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK_CG.txt" if($opt_c);
#	if (not $opt_f) {
		open (my $FINALPOS, ">", $finalPos) or die "Could not open $finalPos: $!";
		open (my $FINALNEG, ">", $finalNeg) or die "Could not open $finalNeg: $!";
		open (my $FINALPOS_NOPEAK, ">", "$finalPos_NOPEAK") or die "Could not open $finalPos: $!";
		open (my $FINALNEG_NOPEAK, ">", "$finalNeg_NOPEAK") or die "Could not open $finalNeg: $!";
		for(my $strand=0; $strand<2; $strand++)
		{
			my $filePos = defined($opt_c) ? "$outDir/$gene\_CG_POS" : "$outDir/$gene\_POS";
			my $fileNeg = defined($opt_c) ? "$outDir/$gene\_CG_NEG" : "$outDir/$gene\_NEG";
			my $fileLast = $strand == 0 ? "$filePos.tsv" : "$fileNeg.tsv";
			LOG($outLog, "\t- Doing gene $CY$gene$N ($fileLast)\n");
			LOG($outLog, "$fileLast doesn't exist!\n") and next if not -e $fileLast;
			LOG($outLog, "$fileLast has no line!\n") and next if `wc -l < $fileLast` == 0;
			open (my $fileLastIn, "<", $fileLast) or LOG($outLog, "Cannot open $fileLast\n" and die "Cannot open $fileLast: $!\n");
			my ($fileLastCount) = `wc -l $fileLast` =~ /^(\d+) /;
			
			my %peak;
			my $start = 1;
			my $end = $minPeakLength;
			my $LCOUNT = 0;
			my $LENGTH = 0; my $NAME = "NA";
			while (my $lineFinal = <$fileLastIn>)
			{
				$LCOUNT ++;
				chomp($lineFinal);
				LOG($outLog, "\tDone $GN$LCOUNT$N\n") if $LCOUNT % 500 eq 0;
				print "\t$fileLast: Done $GN$LCOUNT$N / $fileLastCount\n" if $LCOUNT % 500 eq 0;
				next if $lineFinal =~ /^\#READ/;
				my ($name, @fields) = split("\t", $lineFinal);
				my $readLength = @fields;
				if ($readLength != $LENGTH and $NAME ne "NA") {
					print "PREVIOUS: name=$CY$NAME$N, length=$LGN$LENGTH$N\n";
					print "CURRENT : name=$CY$name$N, length=$LGN$readLength$N\n";
				}
				$LENGTH = $readLength; $NAME = $name;
	
				$peak{$name}{length} = 0;
				#0= not converted (grey) # same
				#1= converted (green)   # same
				#2= non-C (white)       # same
				#3= non-converted CpG (blue) # is non converted
				#4= converted CpG  (black) # is converted
				#5= converted CpG in R-loop (purple) # conv CPG
				#6= no data
				#9= converted C in R-loop (red) # not exist
				# conversion of mine into Jenna's: 3 (no data) becomes 2 (non C)
				
				# peak call based on threshold
				my @peak; my $ispeak = 0;
				for (my $i = 0; $i < @fields - $minPeakLength; $i++)
				{
				#	print "8. GENE=$gene: name=$name, Doing Peak Call!\n" if $name eq "SEQ_119704"; # if $i > 1400 and $i < 1600;
					$peak{$name}{length} ++ if $fields[$i] =~ /^[1459]$/;
					my $chunk = join("", @fields[$i..($i+$minPeakLength)]);
					# calculate % conv from start to end, where end is minPeakLength (e.g. 100 then 1 to 100)
					
					# calculate conC (converted C) and nonConC (non converted)
					# calculate % conv
					# if too few C, (Alex Meissner uses 5 as threshold) we can't call it as anything ... Meissner et al (Nature 2012)
					# if % conv is more than thershold, then turn 1 (conv C) into 9 and 4 (conv C in CpG) into 5
					my ($conC)    = $chunk =~ tr/1459/1459/;
					my ($nonConC) = $chunk =~ tr/03/03/;
					
					my $conPer = $conC + $nonConC >= 5 ? $conC / ($conC + $nonConC) : 0;
					if($conPer >= $opt_t)
					{
						my $checkPeak = 0;
#= CUT1
#= CUT1
						if ($checkPeak == 0) {
							$ispeak = 1;
							$peak[$i] = 1; # from i to i + minPeakLength is a peak
							for(my $j = $i; $j < $i + $minPeakLength; $j++)
							{
								$fields[$j] = 9 if $fields[$j] == 1;
								$fields[$j] = 5 if $fields[$j] == 4;
							}
						}
						else {
							$peak[$i] = 0;
							$ispeak = 0;
						}
					}
					else {
						$peak[$i] = 0;
					}
				}
				$peak{$name}{field} = join("\t", @fields);
				
				# if it's a read that has peak, then we need to put the peak length etc so it'll be presorted for R
				if ($ispeak == 1) {
					my ($maxLength) = 0; my $length = 0; my $lastPos = -2;
					for (my $i = 0; $i < @peak; $i++) {
						if ($lastPos != -2 and $lastPos == $i - 1 and $peak[$i] == 1) {
							$length ++; $lastPos = $i;
							$maxLength = $length if $maxLength < $length;
						}
						else {
							$length = 0; $lastPos = -2;
						}
					}
					$peak{$name}{length} = $maxLength;
					$peak{$name}{peak} = 1;
				}
				else {$peak{$name}{peak} = 0;}
			}
			my $peakcount = 0; my $nopeak = 0; my $total = 0; my $maxAllLength = 0;
			foreach my $name (sort {$peak{$b}{peak} <=> $peak{$a}{peak} || $peak{$b}{length} <=> $peak{$a}{length}} keys %peak)
			{
				#last if $peak{$name}{peak} == 0 and $nopeak > 1000;
				#$peak{$name}{field} =~ tr/19/91/;
				print $FINALPOS "$name\t$peak{$name}{field}\n"        if $strand == 0 and $peak{$name}{peak} == 1;
				print $FINALPOS_NOPEAK "$name\t$peak{$name}{field}\n" if $strand == 0 and $peak{$name}{peak} == 0;
				print $FINALNEG "$name\t$peak{$name}{field}\n"        if $strand == 1 and $peak{$name}{peak} == 1;
				print $FINALNEG_NOPEAK "$name\t$peak{$name}{field}\n" if $strand == 1 and $peak{$name}{peak} == 0;
				$peakcount ++ if $peak{$name}{peak} != 0;
				$nopeak ++ if $peak{$name}{peak} == 0;
				$total ++;
				$maxAllLength = $peak{$name}{length} if $maxAllLength < $peak{$name}{length};
			}
			LOG($outLog, "$finalPos\tlength=$maxAllLength\tpeak=$peakcount\tnopeak=$nopeak\ttotal=$total\n") if $strand == 0;
			LOG($outLog, "$finalNeg\tlength=$maxAllLength\tpeak=$peakcount\tnopeak=$nopeak\ttotal=$total\n") if $strand == 1;
		}
		close $FINALPOS; close $FINALNEG;
		close $FINALPOS_NOPEAK;
		close $FINALNEG_NOPEAK;
#	}	
	# Do edge clearing

PART6END
-e -f -B -D


if ($opt_e) {
	#0 = white
	#1 = intron line
	#2 = exon line
	LOG($outLog, "\t${BR}2a. Parsing in exon intron data from $CY$exonFile$N:\n");
	foreach my $gene (sort keys %{$SEQ}) {
		next if -e "$outFolder/exon/$gene.exon"; # DELETE THIS
		my ($genepos) = `grep -i -P \"\\t$gene\\t\" $geneIndexFile`; chomp($genepos);
		my $length_seq = $SEQ->{$gene}{geneL};
		LOG($outLog, "\n$LRD!!!$N\tWARNING: No sequence for gene $CY$gene$N is found in -s $CY$seqFile$N!\n\n") if $length_seq == 0;
		myFootLib::parseExon($exonFile, $genepos, $gene, $outFolder, $length_seq);
		LOG($outLog, "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n");
	}
}

sub convert_seq {
	my @seq = @{$_[0]};
	my %new;
	for (my $strand = 0; $strand < 2; $strand++) {
		my @new;
		for (my $i = 0; $i < @seq; $i++) {
			my $lastpos = $strand == 0 ? 0 : @seq-1;
			my $add     = $strand == 0 ? 1 : -1; 
			my $nuc1 = $strand == 0 ? "C" : "G";
			my $nuc2 = $strand == 0 ? "G" : "C";
			
			if ($i == 0 || $i == @seq - 1) {
				$new[$i] = "B";
			}
			elsif ($seq[$i] ne $nuc1) {
				$new[$i] = "N";
			}
			else {
				if ($seq[$i+$add] eq $nuc2) {
					$new[$i] = "P" if defined($opt_c);
					$new[$i] = "N" if not defined($opt_c);
				}
				else {
					$new[$i] = "C";
				}
			}
		}
		@{$new{$strand}} = @new;
	}
	return(\%new);
}

=BOX
	my %box; 
	if (defined $opt_B) {
		my $boxFile = $opt_B;
		die "\n${LRD}ERROR!$N: $LCY-B$N <boxfile.bed> does not exist!\n" if not -e $boxFile;
		my @lines = `cat $boxFile`;
		die "\n${LRD}ERROR!$N: $LCY-B$N <boxfile.bed> is empty!!\n" if @lines == 0;
		my $lineCount = 0;
		foreach my $line (@lines) {
			chomp($line);
			next if $line =~ /^\#/;
			my ($chr, $beg, $end, $gene, $val, $strand, @others) = split("\t", $line);
			$lineCount ++;
			$gene = "GENE_$lineCount" if not defined $gene;
			$strand = "+" if not defined $strand;
			$val = 0 if not defined $val;
			my $others = @others == 0 ? "" : "\t" . join("\t", @others);
			$chr=uc($chr);
			push(@{$box{$chr}}, "$chr\t$beg\t$end\t$gene\t$val\t$strand$others");
		}
	}
	return(\%box);


__END__
		if ($line =~ /^Total methylated C's in CpG context:/) {
			my ($num) = $line =~ /Total methylated C's in CpG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CpG from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "meth_CpG");
		}
		if ($line =~ /^Total methylated C's in CHG context:/) {
			my ($num) = $line =~ /Total methylated C's in CHG context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHG from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "meth_CHG");
		}
		if ($line =~ /^Total methylated C's in CHH context:/) {
			my ($num) = $line =~ /Total methylated C's in CHH context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in CHH from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "meth_CHH");
		}
		if ($line =~ /^Total methylated C's in Unknown context:/) {
			my ($num) = $line =~ /Total methylated C's in Unknown context:[ \t]+(\d+\.?\d*)\%[ \t]*$/;
			LOG($outLog, "footLoop.pl::parse_bismark_report failed to get Total methylated C in Unknown from log $LCY$bismark_report_file$N!\n") if not defined $num;
			push(@report, $num);
			push(@header, "meth_Unknown");
		}

