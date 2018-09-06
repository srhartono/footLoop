#!/usr/bin/perl
# V3.7c

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_w $opt_g $opt_L $opt_G $opt_v $opt_n $opt_r $opt_R $opt_B $opt_c); #v $opt_x $opt_R $opt_c $opt_t $opt_n);
getopts("n:vg:w:G:r:R:B:cL:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
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

die "\nUsage: $YW$0$N -n$LCY <footPeak output directory>$N

${LGN}Optionals$N:

-G $LGN<gene to process>$N] 

-c: This is CpG (in vitro)


-r Option to run R scripts $LCY(PNG)$N
   -r 0: do not run any R scripts
   -r 1: run only relevant R scripts (default)
   -r 2: run ALL R scripts
	-r 3: PEAK
	-r 4: NOPK
	-r 5: PEAK_TEMP
	-r 6: NOPK_TEMP

-R Option to run R scripts $LCY(PDF)$N
   -R 0: do not run any R scripts (default)
   -R 1: run only relevant R scripts
   -R 2: run ALL R scripts
	-R 3: PEAK
	-R 4: NOPK
	-R 5: PEAK_TEMP
	-R 6: NOPK_TEMP

-L: subsampling number of lines

-B: <BED3 file> add option to add box in the graph

" if not defined $opt_n;
die "\nERROR: -n footPeak dir $LCY$opt_n$N doesn't exists!\n\nUsage: $YW$0$N -n <footPeak output directory>\n\n" if not -d $opt_n;
die "\nERROR: -L <subsample number>\n" if not defined $opt_L or (defined $opt_L and $opt_L !~ /^\d+$/);
my $subsample = $opt_L;

my ($currMainFolder) = `pwd`; chomp($currMainFolder);
$opt_r = 1 if not defined $opt_r;
$opt_R = 0 if not defined $opt_R;
#my $toggleCpG = defined $opt_c ? 1 : 0;

die "\nERROR: -r has to be 0, 1, or 2! (current: $opt_r)\n\n" if defined $opt_r and $opt_r !~ /^[0123456]$/;
die "\nERROR: -R has to be 0, 1, or 2! (current: $opt_R)\n\n" if defined $opt_R and $opt_R !~ /^[0123456]$/;

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
	my ($resDir) = @_;
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
	open (my $outLog, ">", "$resDir/footPeak_graph_subsample_logFile.txt") or die "\n\nFailed to write to $resDir/footPeak_graph_subsample_logFile.txt: $!\n\n";

	############
	# LOG START
	LOG($outLog, ">footPeak_graph.pl version $version\n");
	LOG($outLog, ">UUID: $uuid\n", "NA");
	LOG($outLog, ">Date: $date\n", "NA");
	LOG($outLog, ">Run script: $0 -n $opt_n -L $opt_L\n", "NA");
	
	############

	my %coor;
	my @lines   = `cat $footPeak_logFile`;
	LOG($outLog, " Parsing footpeak logfile $footPeak_logFile\n");
	my ($label) = `cat $resDir/.LABEL`; chomp($label);
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_logFile!\n\n") if not -e $footPeak_logFile;
	DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $resDir/.LABEL!\n\n") if not -e "$resDir/.LABEL";
	my ($thres, $window);
	foreach my $line (@lines) {
		chomp($line);
		if ($line =~ /^[ \t]*def=.+, coor=.+/) {
			$line =~ s/(^\s+|\s+$)//g;
			my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
			LOG($outLog, "gene=$gene,chr=$CHR,beg=$BEG,end=$END,gene=$GENE,val=$VAL,strand=$STRAND\n");
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
	my @types = qw(CH CG GH GC);
	foreach my $GENE (sort keys %coor) {
		my $mygene = $GENE;
		my @strands = qw(Pos Neg Unk);
		#my $strand = $coor{$GENE}{STRAND};
		for (my $h1 = 0; $h1 < @strands; $h1++) {
			for (my $h2 = 0; $h2 < 4; $h2++) {
				my $strand = $strands[$h1];
				my $type = $types[$h2];
				my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK";
				LOG($outLog, "label=$label, gene$mygene, strand=$strand, peak=$peakFile\n");
				$files{$peakFile} = $mygene;
			}
		}
	}
	my %Rscripts; 
	my $lastfile = -1; #debug
	DIELOG($outLog, "\n\nERROR: There is no files in $LCY\%files$N defined!\n") if (keys %files) == 0;
	my $fileCount = 0;
	my $totalFile = (keys %files);
	my $lastGENE = -1;

	my $currfileCount = 0;
	foreach my $file (sort keys %files) {
		$fileCount ++;
		my $GENE = $files{$file};
		my $STRAND = $coor{$GENE}{STRAND};
		my $geneStrandPrint = $STRAND eq "Pos" ? "$LRD$STRAND$N" : $STRAND eq "Neg" ? "$LCY$STRAND$N" : "$LGN$STRAND$N";

   	if (defined $opt_G and $file !~ /$opt_G/i) {
   	   LOG($outLog, date() . " Skipped $LCY$file$N as it doesn't contain $LGN-G $opt_G$N\n");
   	   next;
   	}
		if ($GENE ne $lastGENE) {
			LOG($outLog, "\n$YW -------- $fileCount/$totalFile Doing$LPR $GENE$N (STRAND = $geneStrandPrint) ---------$N\n\n") if $GENE ne $lastGENE;
			$currfileCount = 0;
		}
		$lastGENE = $GENE;
		if (defined $opt_g and $GENE ne $opt_g) {
			LOG($outLog, date() . " $LCY Skipped $GENE$N (requested gene is $LGN$opt_g$N\n");
			next;
		} 
		if (defined $opt_w and $file !~ /$opt_w/) {
			LOG($outLog, date() . " Skipped $LCY$file$N (requested want is $LGN$opt_w$N\n");
			next;
		} 

		$lastfile = $file;
		next if not defined $files{$file};

		my ($pk_filename) =  getFilename($file, 'full') . ".out";
		my $peakFile      =  "$resDir/.CALL/$pk_filename";
		my $nopkFile      =  $peakFile;
		my $cluster_file  = "$resDir/FOOTCLUST/CLUST_LOCAL/" . getFilename($file, "full") . ".local.bed.indiv.clust";
		my $kmer_file     =  "$resDir/FOOTCLUST/.TEMP/$pk_filename";
			$nopkFile      =~ s/\.PEAK.out$/.NOPK.out/;
			$kmer_file     =~ s/.out$/.local.bed.clust.kmer/;
		my $bedFile       =  "$resDir/PEAKS_LOCAL/$pk_filename.local.bed";
			$bedFile       =~ s/.out.local.bed/.local.bed/;
		my $totpeak = -e $peakFile ? linecount($peakFile) : 0;
		my $totnopk = -e $nopkFile ? linecount($nopkFile) : 0;
		my ($type) = $peakFile =~ /_(CH|CG|GH|GC)./;
		my $parseName = parseName($pk_filename);
		my ($label2, $gene2, $strand2, $window2, $thres2, $type2) = @{$parseName->{array}};
		DIELOG($outLog, date() . "Undefined label2=$label2, gene2=$gene2, strand2=$strand2, window2=$window2, thres2=$thres2, type2=$type2 from parsing filename=$LCY$pk_filename$N\n") if not defined $label2 or not defined $gene2 or not defined $strand2 or not defined $window2 or not defined $thres2 or not defined $type2;

		for (my $p = 0; $p < 2; $p ++) {
			$currfileCount ++;
			my $currFile = $p == 0 ? $peakFile : $nopkFile;
			my $currFileID = $currFile . ".id";
			my $curr_cluster_file = $cluster_file;
				$curr_cluster_file =~ s/PEAK/NOPK/ if $p != 0;
			my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
			my ($geneStrand, $readStrand, $rconvType) = ($STRAND, $strand2, $type2);
			my $readStrandPrint = $readStrand eq "Pos" ? "$LRD$readStrand$N" : $readStrand eq "Neg" ? "$BU$readStrand$N" : "$LPR$readStrand$N";
			my $rconvTypePrint  = $rconvType =~ /^(CG|GC)$/ ? "$YW$rconvType$N" : $rconvType;
			my $max_cluster = 0;
			if (-e $curr_cluster_file and -s $curr_cluster_file > 0) {
				my @curr_cluster = `cut -f6 $curr_cluster_file|grep -v clust`; chomp(@curr_cluster);
				foreach my $curr_cluster (sort {$b <=> $a} @curr_cluster) {
					$max_cluster = $curr_cluster; last;
				}
			}
			$max_cluster = 0 if not defined $max_cluster;
			my $summary = $p == 0 ? "PEAK: $LGN$totpeak$N, cluster=$LCY$max_cluster$N" : "NOPK: $LGN$totnopk$N";
			#$currFilename =~ s/\.0_orig_//i;
			#my $flag = $readStrand eq "Unk" ? "" : $p == 0 ? $PEAK->[0] : $PEAK->[1]; #PEAK : NOPK;
			my $flag = getFlag($currFile, $geneStrand, $readStrand, $rconvType);
			my $pngoutDir = $flag;
			my $pdfoutDir = $flag;
			LOG($outLog, date() . "$LGN$currfileCount.$N $flag $readStrandPrint $rconvTypePrint $LCY$currFile$N\n");
			LOG($outLog, "\t\tCurrfile           = $LCY$currFile$N
\t\tcurr_cluster_file  = $LCY$curr_cluster_file$N
\t\tkmer_File          = $LPR$kmer_file$N
\t\tbedFile            = $BU$bedFile$N
\t\ttotpeak = $LGN$totpeak$N, nopk = $LGN$totnopk$N
","NA");
			LOG($outLog, "$readStrandPrint\t$rconvTypePrint\t$flag\n","NA");

####SPACESHIP
			my $c = defined $opt_c ? "_C" : "";
			my $resDir2 = getFullpath($resDir);
			my $madePNG = 0; my $madePDF = 0;
			$madePNG = 1 if (  
					($opt_r eq 1 and $flag =~ /^(PEAK$c|NOPK$c|PEAK_TEMP$c|NOPK_TEMP$c)$/) or 
					($opt_r eq 2) or
					($opt_r eq 3 and $flag eq "PEAK$c") or
					($opt_r eq 4 and $flag eq "NOPK$c") or
					($opt_r eq 5 and $flag eq "PEAK_TEMP$c") or
					($opt_r eq 6 and $flag eq "NOPK_TEMP$c")
			);
			$madePDF = 1 if (  
					($opt_R eq 1 and $flag =~ /^(PEAK$c|NOPK$c|PEAK_TEMP$c|NOPK_TEMP$c)$/) or 
					($opt_R eq 2) or
					($opt_R eq 3 and $flag eq "PEAK$c") or
					($opt_R eq 4 and $flag eq "NOPK$c") or
					($opt_R eq 5 and $flag eq "PEAK_TEMP$c") or
					($opt_R eq 6 and $flag eq "NOPK_TEMP$c")
			);
			LOG($outLog, date() . " --> DEBUG flag=$flag madePNG = $madePNG\n","NA");
#			my $madePDF = ($opt_R == 0) ? 0 : ($opt_R == 1 and defined $opt_c and $flag =~ /^(PEAK_C|NOPK_C|PEAK_TEMP_C|NOPK_TEMP_C)$/) ? 1 : ($opt_R == 1 and $flag =~ /(ALL|RCONV|_C)/) ? 0 : 1;
			LOG($outLog, date() . " --> DEBUG flag=$flag madePDF = $madePDF\n","NA");
			my $pngout = "$resDir/PNG/$pngoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.png";
			my $pdfout = "$resDir/PDF/$pdfoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.pdf";
			my $lenpdfout = "$resDir/PDF/$pngoutDir/ALL/$currFilename\_length.$flag.rand.$subsample.pdf";

			$scp{"scp $user\@crick.cse.ucdavis.edu:$resDir2/PNG/$pngoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.png ./"} = 1 if $madePNG eq 1;
			$scp{"scp $user\@crick.cse.ucdavis.edu:$resDir2/PNG/$pngoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.png.c_conv.png ./"} = 1 if $madePNG eq 1;
			$scp{"scp $user\@crick.cse.ucdavis.edu:$resDir2/PDF/$pdfoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.pdf ./"} = 1 if $madePDF eq 1;
			$scp{"scp $user\@crick.cse.ucdavis.edu:$resDir2/PDF/$pdfoutDir/ALL/$currFilename.$flag.rand.$subsample.heatmap.pdf.c_conv.pdf ./"} = 1 if $madePDF eq 1;

#			$scp{"scp $user\@crick.cse.ucdavis.edu:$resDir2/PDF/$pdfoutDir/ALL/$currFilename.$flag.rand.$subsample.pdf ./"} = 1 if $madePDF eq 1;
			#LOG($outLog, "!HERE STRAND=$STRAND strand=$strand type=$type label=$label pngoutDir = $pngoutDir, p = $p, PNGOUT = $pngout\n");
			my ($RscriptPDF, $RscriptPNG, $RscriptPDF_nopk_ALL, $RscriptPNG_nopk_ALL);
			my $Rscript = "
.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\", \"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\", .libPaths()) )
library(labeling)\nlibrary(ggplot2)\nlibrary(reshape2)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(RColorBrewer)
";
			
			#if (-e $currFile and linecount($currFile) > 10) {
			my $totread = $totpeak + $totnopk;
			if (not -e $currFile or (-e $currFile and linecount($currFile) < 3)) {
				$Rscript .= "MYPNGDESC \nplot(NA,xlim=c(1,100),ylim=c(1,100),xaxt=F,xlab=NA,ylab=NA,bty=\"n\")\ntext(50,50,cex=3,labels=c(\"$currFilename\n\nPEAK = $totpeak / $totread\"))\ndev.off()\n";
				$RscriptPNG = $Rscript; $RscriptPNG =~ s/MYPNGDESC/png\("$pngout",1000,1000\)/g;
				$RscriptPDF = $Rscript; $RscriptPDF =~ s/MYPNGDESC/pdf\("$pdfout",height=7,width=7\)/g;
				$RscriptPNG_nopk_ALL = $RscriptPNG;
				$RscriptPDF_nopk_ALL = $RscriptPDF;
			}
			else {
				open (my $incurrFile, "<", $currFile) or DIELOG($outLog, "Failed to open $currFile: $!\n");
				open (my $outcurrFileID, ">", "$currFile.id") or DIELOG($outLog, "Failed to open $currFile.id: $!\n");
				while (my $line = <$incurrFile>) {
					chomp($line);
					my ($id) = split("\t", $line);
					my ($num1, $num2, $num3) = $id =~ /^.*m(\d+_\d+)_.+\/(\d+)\/(ccs|\d+_\d+)/;
					   ($num1, $num2, $num3) = $id =~ /^.*m(\d+_\d+)\/(\d+)\/(ccs|\d+_\d+)/ if not defined $num1;
					DIELOG($outLog, "\n\n" . date() . " Failed to parse numbers from currFile=$LCY$currFile$N id = $LPR$id$N\n\n") if not defined $num1 or not defined $num2 or not defined $num3;
					$num3 = 0 if $num3 eq "ccs";
					my $num = "$num1$num2$num3";
			      $num =~ s/_//g;
					print $outcurrFileID "$id\t$num\n";
				}
				close $incurrFile;
				close $outcurrFileID;
				my ($R) = Rscript($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout, $resDir, $subsample);

				# Read Table
				$Rscript .= $R->{readTable};
				if (-e $curr_cluster_file and -s $curr_cluster_file > 0 and $currFile !~ /\.NOPK\./) {
					LOG($outLog, "\n\n$LRD CLUSTER $N\n\n");
					$Rscript .= $R->{clusterFile};
				} 
				else {
					$Rscript .= $R->{noclusterFile};
				}
				$Rscript .= $R->{mainplot_nopk_rand_subsample};
				if (-e $bedFile and linecount($bedFile) > 0 and $currFile !~ /\.NOPK\./) {
					$Rscript .= $R->{peakbedFile};
				}
				$Rscript .= $R->{secondplotConversionGraph_rand_subsample};
				$Rscript .= $R->{Scale};
				$RscriptPNG .= $Rscript . $R->{PNG_nopk_rand_subsample};
				$RscriptPDF .= $Rscript . $R->{PDF_nopk_rand_subsample};
			}

			open (my $outRscriptPNG, ">", "$currFile.PNG.$subsample.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PNG.$subsample.R: $!\n") and print $outLog $Rscript and next);
			print $outRscriptPNG $RscriptPNG;

			$Rscripts{"$currFile.PNG.$subsample.R"}{summary} = $summary;
			$Rscripts{"$currFile.PNG.$subsample.R"}{runR} = (defined $opt_c and $flag =~ /^(NOPK_C|PEAK_C|NOPK_TEMP_C|PEAK_TEMP_C)$/) ? 1 : $flag =~ /(ALL|RCONV|_C)/ ? 0 : 1;
			$Rscripts{"$currFile.PNG.$subsample.R"}{runPNG} = $madePNG;
			$Rscripts{"$currFile.PNG.$subsample.R"}{runPDF} = $madePDF;
			$Rscripts{"$currFile.PNG.$subsample.R"}{runType} = $flag;
			close $outRscriptPNG;

			open (my $outRscriptPDF, ">", "$currFile.PDF.$subsample.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PDF.$subsample.R: $!\n") and print $outLog $Rscript and next);
			print $outRscriptPDF $RscriptPDF;
			close $outRscriptPDF;

#			if ($currFile =~ /\.NOPK\./) {
#				open (my $outRscriptPNG_nopk_ALL, ">", "$currFile.PNG_nopk_ALL.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PNG_nopk_ALL.R: $!\n") and print $outLog $Rscript and next);
#				print $outRscriptPNG_nopk_ALL $RscriptPNG_nopk_ALL;
#				$Rscripts{"$currFile.PNG_nopk_ALL.R"}{summary} = $summary;
#				$Rscripts{"$currFile.PNG_nopk_ALL.R"}{runR} = 0;#$flag =~ /(ALL|RCONV|_C)/ ? 0 : 1;
#				$Rscripts{"$currFile.PNG_nopk_ALL.R"}{runType} = $flag . "_ALL";
#				close $outRscriptPNG_nopk_ALL;
#
#				open (my $outRscriptPDF_nopk_ALL, ">", "$currFile.PDF_nopk_ALL.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PDF_nopk_ALL.R: $!\n") and print $outLog $Rscript and next);
#				print $outRscriptPDF_nopk_ALL $RscriptPDF_nopk_ALL;
#				close $outRscriptPDF_nopk_ALL;
#			}
		}
	}
	LOG($outLog, "\n\n$YW ----------------- Running R Scripts ------------------$N\n\n");
	$fileCount = 0;
	$totalFile = (keys %Rscripts);
	
	# open outRscripts for Rscripts that aren't relevant
	open (my $outR_notrelevantPNG, ">", "$resDir/PNG/footPeak_graph_Rscripts.PNG.$subsample.sh") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/PNG/footPeak_graph_Rscripts.PNG.$subsample.sh: $!\n");
	open (my $outR_notrelevantPDF, ">", "$resDir/PDF/footPeak_graph_Rscripts.PDF.$subsample.sh") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/PDF/footPeak_graph_Rscripts.PDF.$subsample.sh: $!\n");
	foreach my $outRscriptPNG (sort keys %Rscripts) {
		my $summary = $Rscripts{$outRscriptPNG}{summary};
		my $runR = $Rscripts{$outRscriptPNG}{runR};
		my $runPNG = $Rscripts{$outRscriptPNG}{runPNG};
		my $runPDF = $Rscripts{$outRscriptPNG}{runPDF};
		my $runType = $Rscripts{$outRscriptPNG}{runType};
		my $outRscriptPDF = $outRscriptPNG; 
			$outRscriptPDF =~ s/PNG.$subsample.R$/PDF.$subsample.R/g;
			$outRscriptPDF =~ s/PNG_nopk_ALL.$subsample.R$/PDF_nopk_ALL.$subsample.R/g;
		LOG($outLog, date() . " $LCY Skipped $outRscriptPNG$N (requested gene is $LGN$opt_g$N\n") and next if defined $opt_g and $outRscriptPNG !~ /$opt_g/;
		$fileCount ++;
		my $RLOG = 0;
		
		#if (($opt_r == 2) or ($opt_r == 1 and $runR == 1)) {
		if ($runPNG eq 1) {
			LOG($outLog, "\n" . date() . "flag=$LPR$runType$N, $LGN$fileCount/$totalFile$N. Running$LGN PNG$N $LCY$outRscriptPNG$N: $summary\n");
			LOG($outLog, date() . "\tprinted to ${LCY}footPeak_graph_Rscripts.sh$N: cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPNG > $outRscriptPNG.LOG 2>&1\n","NA");
			print $outR_notrelevantPNG "cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPNG > $outRscriptPNG.LOG 2>&1 #$runType\n";
			$RLOG = system("cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPNG > $outRscriptPNG.LOG 2>&1");
		}
		else {
			LOG($outLog, "\n" . date() . "flag=$LPR$runType$N, $LGN$fileCount/$totalFile$N.$LRD Not$N running$LGN PNG$N $LCY$outRscriptPNG$N: $summary\n");
			LOG($outLog, date() . "\tprinted to ${LCY}footPeak_graph_Rscripts.sh$N: cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPNG > $outRscriptPNG.LOG 2>&1\n","NA");
			print $outR_notrelevantPNG "cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPNG > $outRscriptPNG.LOG 2>&1 #$runType\n";
			$RLOG = 1;
		}
		my $prevRLOG = $RLOG;
#		if (($opt_r == 2) or ($opt_r == 1 and $runR == 1)) {
		if ($runPNG eq 1) {
			if ($RLOG ne 0) {
				if (not -e "$outRscriptPNG.LOG") {
					$RLOG = "\t$outRscriptPNG.LOG cannot be found!\n" if $runR == 1;
				}
				else {
					my @RLOG = `tail -n 5 $outRscriptPNG.LOG`;
					$RLOG = "\t" . join("\n\t", @RLOG);
				}
				LOG($outLog, date() . "\t--> ${LRD}Failed$N to run_Rscript.pl $outRscriptPNG: $prevRLOG, LOG:\n$RLOG\n") if $runR == 1;
			}
			else {
				LOG($outLog, date() . "\t--> ${LGN}Success$N on running run_Rscript.pl $LCY$outRscriptPNG$N\n");
			}
		}

		$RLOG = 0;
		#if (($opt_R == 2) or ($opt_R == 1 and $runR == 1)) {
		if ($runPDF eq 1) {
			LOG($outLog, "\n" . date() . "flag=$LPR$runType$N, $LGN$fileCount/$totalFile$N. Running$LGN PDF: $LCY$outRscriptPDF$N: $summary\n");
			LOG($outLog, date() . "\tprinted to ${LCY}footPeak_graph_Rscripts.sh$N: cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPDF > $outRscriptPDF.LOG 2>&1\n","NA");
			print $outR_notrelevantPDF "cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPDF > $outRscriptPDF.LOG 2>&1 #$runType\n";
			$RLOG = system("cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPDF > $outRscriptPDF.LOG 2>&1");
		}
		else {
			LOG($outLog, "\n" . date() . "flag=$LPR$runType$N, $LGN$fileCount/$totalFile$N.$LRD Not$N running$LGN PDF$N $LCY$outRscriptPDF$N: $summary\n");
			LOG($outLog, date() . "\tprinted to ${LCY}footPeak_graph_Rscripts.sh$N: cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPDF > $outRscriptPDF.LOG 2>&1\n","NA");
			print $outR_notrelevantPDF "cd $currMainFolder/ && R --vanilla --no-save < $outRscriptPDF > $outRscriptPDF.LOG 2>&1 #$runType\n";
			$RLOG = 1;
		}
		$prevRLOG = $RLOG;
		#if (($opt_R == 2) or ($opt_R == 1 and $runR == 1)) {
		if ($runPDF eq 1) {
			if ($RLOG ne 0) {
				if (not -e "$outRscriptPDF.LOG") {
					$RLOG = "\t$outRscriptPDF.LOG cannot be found!\n" if $runR == 1;
				}
				else {
					my @RLOG = `tail -n 5 $outRscriptPDF.LOG`;
					$RLOG = "\t" . join("\n\t", @RLOG);
				}
				LOG($outLog, date() . "\t--> ${LRD}Failed$N to run_Rscript.pl $outRscriptPDF: $prevRLOG, LOG:\n$RLOG\n") if $runR == 1;
			}
			else {
				LOG($outLog, date() . "\t--> ${LGN}Success$N on running run_Rscript.pl $LCY$outRscriptPDF$N\n");
			}
		}
	}
#	LOG($outLog, date . "\tcd $resDir && run_Rscript.pl *MakeHeatmap.R\n");
#	system("cd $resDir && run_Rscript.pl *MakeHeatmap.R") if not defined $opt_x and defined $$opt_R;
	LOG($outLog, "\n\n$YW ----------------- SCP PATHS ------------------$N\n\n");
	foreach my $file (sort keys %scp) {
		LOG($outLog, "$file\n");
	}
}

###############
# Subroutines #
###############


sub Rscript {
	my ($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout, $resDir, $subsample) = @_;
	my $currFileID = $currFile . ".id";
	my $R;
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
	($bedFilename) =~ s/(CG|CH|GC|GH).+$/$1/;
	my $peakminVal = $currFile =~ /PEAK.out$/ ? 8 : 6;
	my $plotminReads = 5;

# -------------------- $R->{readTable}
	$R->{readTable} = "	
mywindow = 1000
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
theme_blank\$panel.border= element_blank()#structure(c(0, 0, 0,0), unit = \"lines\", valid.unit = 3L, class = \"unit\")
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




# -------------------- $R->{mainplot_nopk_rand_subsample}
	$R->{mainplot_nopk_rand_subsample} .= "

##########################
# Main Plot Part 1 nopk subset subsample

# Randomly take reads
df.total = dim(df)[1]
subsample = $subsample;
set.seed(42)
if (df.total > $subsample) {
	df.rand.subsample = df[sample( seq(1,df.total) ,$subsample, replace=F),]
	print(\"randoming $subsample. DF TOTAL = \"); print(df.total)
} else {
	df.rand.subsample = df
	print(\"NOT randoming $subsample. DF TOTAL = \"); print(df.total)
}

# sorting by hclust
#if (dim(df.rand.subsample)[1] < 1000) {
#	h = hclust(dist(df.rand.subsample[,-1]))
#	df.rand.subsample = df.rand.subsample[h\$order,]
#} else if (dim(df.rand.subsample)[2] < 10) {
#	mysum = apply(df.rand.subsample[,-1],1,sum)
#	df.rand.subsample = df.rand.subsample[order(mysum),]
#}
df.rand.subsample = df.rand.subsample[order(df.rand.subsample\$y),]
df.rand.subsample\$y = seq(1,dim(df.rand.subsample)[1])

#write.table(df,file=\"$resDir/.CALL/$currFilename.out.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
write.table(df.rand.subsample,file=\"$currFile.$subsample.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
print(\"Wrote to $currFile.$subsample.rand\")

dm.rand.subsample = melt(df.rand.subsample,id.vars=c(\"V1\",\"y\"))
print(dim(df.rand.subsample))
print(dim(dm.rand.subsample))
dm.rand.subsample\$variable = as.numeric(as.character(dm.rand.subsample\$variable))

p.rand.subsample = ggplot(dm.rand.subsample,aes(variable,y)) +  
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

p.heatmaponly.rand.subsample = ggplot(dm.rand.subsample,aes(variable,y)) +  
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

	p.rand.subsample.pdf = p.rand.subsample + theme(plot.title = element_text(size = 10*p.pdf.scale))
	p.rand.subsample.png = p.rand.subsample + theme(plot.title = element_text(size = 10*p.png.scale))

";

	$R->{peakbedFile} = "

##########################
# Peak BedFile
bed = read.table(\"$bedFile\",sep=\"\\t\")
bed = merge(subset(df.rand.subsample,select=c(\"V1\",\"y\")),bed,by=\"V1\")

if (length(bed) != 0 & dim(bed)[1] > 0) {
	bed\$variable=1
	p.heatmaponly.rand.subsample =	p.heatmaponly.rand.subsample + 
			 	geom_rect(
					data=bed, 
					aes(xmin=V2, xmax=V3, ymin=y-0.5, ymax=y+0.5),
					fill=rgb(0,0,0,alpha=0),
					color=rgb(1,0,0,alpha=0.25),
			 		size=0.5*p.png.scale # SIZE
				)
}
";


	$R->{secondplotConversionGraph_rand_subsample} = "

# Calculate % Conversion
df2.rand.subsample = subset(df.rand.subsample,select=c(-V1,-y));
if (length(df2.rand.subsample[df2.rand.subsample < $peakminVal]) > 0) {
	df2.rand.subsample[df2.rand.subsample < $peakminVal] = 0;
}
if (length(df2.rand.subsample[df2.rand.subsample >= $peakminVal]) > 0) {
	df2.rand.subsample[df2.rand.subsample >= $peakminVal] = 1
}
df2.rand.subsample = data.frame(x=seq(1,dim(df2.rand.subsample)[2]), y=apply(df2.rand.subsample,2,mean))
if (dim(df2.rand.subsample[df2.rand.subsample\$y > 0,])[1] > $plotminReads) {
	df2.rand.subsample = df2.rand.subsample[df2.rand.subsample\$y > 0,]
	df2.rand.subsample\$x = as.numeric(as.character(df2.rand.subsample\$x));
	df2.rand.subsample\$y = as.numeric(as.character(df2.rand.subsample\$y));
	df2.rand.subsample\$x2 = df2.rand.subsample\$x
	df2.rand.subsample\$y2 = df2.rand.subsample\$y
	for (i in 1:(dim(df2.rand.subsample)[1]-10)) {
		a = df2.rand.subsample[df2.rand.subsample\$x >= df2.rand.subsample[i,]\$x & df2.rand.subsample\$x <= df2.rand.subsample[i+10-1,]\$x,]
		if (length(a) != 0 & dim(a)[1] != 0) {
			df2.rand.subsample[i,]\$y2 = mean(a\$y)
			df2.rand.subsample[i,]\$x2 = mean(a\$x)
		}
	}
	mins = seq(1,as.integer(df2.rand.subsample[1,]\$x2)-1,10)
	maxs = seq(max(df2.rand.subsample\$x2),dim(df.rand.subsample)[2]-2,10)
	df2.rand.subsample = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2.rand.subsample)
	df2.rand.subsample = rbind(df2.rand.subsample,data.frame(x=maxs,y=0,x2=maxs,y2=0))
} else {
	df2.rand.subsample = data.frame(x=seq(1,dim(df.rand.subsample)[2]), y=0, x2=seq(1,dim(df.rand.subsample)[2]), y2=0);
}
 
# P2 % Conversion XY Plot
p2.png.scale = p.png.scale
p2.pdf.scale = p.pdf.scale

p2.rand.subsample.png = 
	ggplot(df2.rand.subsample,aes(x2,y2)) + 
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

p2.rand.subsample.pdf = 
	ggplot(df2.rand.subsample,aes(x2,y2)) + 
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







# -------------------- $R->{clusterFile}

	$R->{clusterFile} = "

#####################
# Cluster
clust = read.table(\"$curr_cluster_file\",header=F,sep=\"\\t\",colClasses=c(\"factor\",\"integer\",\"integer\",\"integer\",\"integer\",\"integer\",\"factor\"))
colnames(clust) = c(\"id\",\"x\",\"xmax\",\"y\",\"ymax\",\"clust\",\"id2\")
clust = subset(clust,select=-id2)
print(head(clust))
#clust\$id = gsub(\"^(.+)\\\\.[0-9]+\$\",\"\\\\1\",clust\$id,perl=T)
clust\$y = seq(1,dim(clust)[1])
clust2 = as.data.frame(aggregate(clust\$y,by=list(clust\$clust),min))
clust2\$ymax = aggregate(clust\$y,by=list(clust\$clust),max)\$x
clust2\$xpos0 = aggregate(clust\$x,by=list(clust\$clust),min)\$x
clust2\$xpos1 = aggregate(clust\$xmax,by=list(clust\$clust),max)\$x
colnames(clust2) = c(\"clust\",\"ymin\",\"ymax\",\"xpos0\",\"xpos1\")
clust2\$xmin = 1
clust2\$xmax = 70
clust2\$clust = clust2\$clust + 10
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
	#theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank())


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
	#theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank())

";

# -------------------- $R->{noclusterFile}
	$R->{noclusterFile} = " 

##########################
# Cluster: NO CLUSTER
# sorting by hclust
df=subset(df,select=-id)
if (dim(df)[1] < 1000) {
	h = hclust(dist(df[,-1]))
	df = df[h\$order,]
} else if (dim(df)[2] < 10) {
	mysum = apply(df[,-1],1,sum)
	df = df[order(mysum),]
}
df\$y = seq(1,dim(df)[1])

";

# -------------------- $R->{peakbedFile}

#	$R->{peakbedFile} = "
#
###########################
## Peak BedFile
#bed = read.table(\"$bedFile\",sep=\"\\t\")
#bed = merge(subset(df,select=c(\"V1\",\"y\")),bed,by=\"V1\")
#
#";

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
";

	$R->{box} = "

	box = read.table(\"$boxFile\",sep=\"\\t\")
	p.png = p.png + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.pdf = p.pdf + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	
	";

	$R->{box_nopk} = "

	box = read.table(\"$boxFile\",sep=\"\\t\")
	p.rand.subsample.png  = p.rand.subsample.png  + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.subsample\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.rand.subsample.pdf  = p.rand.subsample.pdf  + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.subsample\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	p.rand.1000.png = p.rand.1000.png + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,color=\"black\",size=0.5*p.png.scale) # SIZE
	p.rand.1000.pdf = p.rand.1000.pdf + geom_rect(data=box,aes(x=0,y=0,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,color=\"black\",size=0.5*p.pdf.scale) # SIZE
	
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

#write.table(df,file=\"$resDir/.CALL/$currFilename.out.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
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
#   print(paste(\"i=\",i,\", beg=\",beg,\", end=\",end,sep=\"\"))
   dm.temp = dm[dm\$y >= beg & dm\$y <= end,]
   dm.temp\$y = dm.temp\$y - beg + 1
   print(head(dm.temp))
#	dm.temp = dm[seq( ((i-1)*mywindow+1), (i*mywindow) ),]
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
}
";


# -------------------- $R->{mainplotClusterAddition}
	$R->{mainplotClusterAddition} = "

# Main Plot Cluster Addition
# geom rect default size = 0.5
# geom text default size = 10
p.png = p + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.png.scale) + # SIZE
	 geom_rect(data=clust2,aes(color=as.factor(clust),x=xpos0,y=ymin,xmin=xpos0,xmax=xpos1,ymin=ymin,ymax=ymax),size=0.5*p.png.scale,fill=rgb(1,1,1,alpha=0),lwd=1*p.png.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=clust-10),hjust=0,size=10*p.png.scale) +
	 theme(plot.title = element_text(size = 10*p.png.scale))

p.pdf = p + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.pdf.scale) + # SIZE
	 geom_rect(data=clust2,aes(color=as.factor(clust),x=xpos0,y=ymin,xmin=xpos0,xmax=xpos1,ymin=ymin,ymax=ymax),size=0.5*p.pdf.scale,fill=rgb(1,1,1,alpha=0),lwd=1*p.pdf.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=clust-10),hjust=0,size=20*p.pdf.scale) +
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
if (length(df2[df2 < $peakminVal]) > 0) {
	df2[df2 < $peakminVal] = 0;
}
if (length(df2[df2 >= $peakminVal]) > 0) {
	df2[df2 >= $peakminVal] = 1
}
df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean))
if (dim(df2[df2\$y > 0,])[1] > $plotminReads) {
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
	p2.png.scale = p.png.scale
	p2.pdf.scale = p.pdf.scale

p2.png = 
	ggplot(df2,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_bw() +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme_blank + #(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
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
	theme_blank +#(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
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
if (dim(df2.rand.1000[df2.rand.1000\$y > 0,])[1] > $plotminReads) {
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
p2.png.scale = p.png.scale
p2.pdf.scale = p.pdf.scale

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
png(\"$pngout\",width=totalwidth,height=totalheight)
if (mynrow == 3) {
	grid.arrange(p.png,p2.png,p3.png,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p.png,p2.png,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.ALL.heatmap.png\"
png(pngout_heatmap_only,width=totalwidth,height=totalheight)
print(p.heatmaponly)
dev.off()

# PNG all Conv
pngout_peak_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.ALL.c_conv.png\"
png(pngout_peak_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2.png)
dev.off()

";

	$R->{PNG_nopk} = "

# PNG
totalheight = (dim(df.rand.1000)[1] + 31.25) * myscale
png(\"$pngout\",width=totalwidth,height=totalheight)
grid.arrange(p.rand.1000.png,p2.rand.1000.png,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df.rand.1000)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.ALL.heatmap.png\"
png(pngout_heatmap_only,width=totalwidth,height=totalheight)
print(p.heatmaponly.rand.1000)
dev.off()

# PNG all Conv
pngout_nopk_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.ALL.c_conv.png\"
png(pngout_nopk_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2.png)
dev.off()

";

	$R->{PNG_nopk_rand_subsample} = "

pngout_nopk_rand_subsample = \"$pngoutFolder/ALL/$pngoutFilename.rand.subsample.png\"
# PNG
#totalheight = (dim(df.rand.subsample)[1] + 31.25) * myscale
#png(pngout_nopk_rand_subsample,width=totalwidth,height=totalheight)
#grid.arrange(p.rand.subsample.png,p2.rand.subsample.png,ncol=1,nrow=mynrow,heights=totalratio)
#dev.off()

# PNG HEATMAP ONLY
totalheight = dim(df.rand.subsample)[1] * myscale
pngout_heatmap_only = \"$pngoutFolder/$pngoutFilename\"
png(pngout_heatmap_only,width=totalwidth,height=totalheight)
print(p.heatmaponly.rand.subsample)
dev.off()

# PNG all Conv
pngout_nopk_all_c_conv = \"$pngoutFolder/$pngoutFilename.c_conv.png\"
png(pngout_nopk_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2.rand.subsample.png)
dev.off()

";

	$R->{PNG_nopk_ALL} = "

#plot_list = list();
for (i in seq(1,as.integer(dim(df)[1] / mywindow) + 1)) {
	pngout_nopk = \"$pngoutFolder/ALL/$pngoutFilename\"
	pngout_nopk = gsub(\"^(.+).png\$\",\"\\\\1\",pngout_nopk,perl=T)
	pngout_nopk = paste(pngout_nopk,\".\",i,\".png\",sep=\"\")
	if (i == max(seq( 1,as.integer(dim(df)[1]/mywindow) + 1 ))) {
		currtotalheight = totalheight_nopk_last
		currtotalratio = totalratio_nopk_last
		print(paste(i,currtotalheight,currtotalratio))
		png(pngout_nopk,width=totalwidth,height=currtotalheight)
		grid.arrange(plot_list[[i]],p2,ncol=1,nrow=mynrow,heights=currtotalratio)
		dev.off()
	} else {
		currtotalheight = totalheight_nopk
		print(paste(i,currtotalheight))
		png(pngout_nopk,width=totalwidth,height=currtotalheight)
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
pdf(\"$pdfout\",width=currwidth,height=currheight)
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
pdf(pdfout_heatmap_only,width=currwidth,height=currheight)
print(p.heatmaponly)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_peak_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_peak_all_c_conv,width=currwidth,height=currheight)
grid.arrange(p2.pdf)
dev.off()

";

	$R->{PDF_nopk} = "

# PDF
currheight = (dim(df.rand.1000)[1] + 31.25) * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdf(\"$pdfout\",width=currwidth,height=currheight)
grid.arrange(p.rand.1000.pdf,p2.rand.1000.pdf,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()


# PDF HEATMAP ONLY
currheight = dim(df.rand.1000)[1] * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_heatmap_only = \"$pdfoutFolder/$pdfoutFilename.ALL.heatmap.pdf\"
pdf(pdfout_heatmap_only,width=currwidth,height=currheight)
print(p.heatmaponly.rand.1000)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_all_c_conv = \"$pdfoutFolder/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=currwidth,height=currheight)
grid.arrange(p2.pdf)
dev.off()

";

	$R->{PDF_nopk_rand_subsample} = "
currheight = (dim(df.rand.subsample)[1] + 31.25) * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_rand_subsample = \"$pdfoutFolder/ALL/$pdfoutFilename.rand.subsample.pdf\"

# PDF
#pdf(pdfout_nopk_rand_subsample,width=currwidth,height=currheight)
#grid.arrange(p.rand.subsample.pdf,p2.rand.subsample.pdf,ncol=1,nrow=mynrow,heights=totalratio)
#dev.off()

# PDF HEATMAP ONLY
currheight = dim(df.rand.subsample)[1] * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_heatmap_only = \"$pdfoutFolder/$pdfoutFilename\"
pdf(pdfout_heatmap_only,width=currwidth,height=currheight)
print(p.heatmaponly.rand.subsample)
dev.off()


# PDF all Conv
currheight = 31.25 * myscale / $PDFSCALE
currwidth = totalwidth / $PDFSCALE
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=currwidth,height=currheight)
grid.arrange(p2.rand.subsample.pdf)
dev.off()

";

	$R->{PDF_nopk_ALL} = "

#plot_list = list();
for (i in seq(1,as.integer(dim(df)[1] / mywindow) + 1)) {
	pdfout_nopk = \"$pdfoutFolder/ALL/$pdfoutFilename\"
	pdfout_nopk = gsub(\"^(.+).pdf\$\",\"\\\\1\",pdfout_nopk,perl=T)
	pdfout_nopk = paste(pdfout_nopk,\".\",i,\".pdf\",sep=\"\")
	if (i == max(seq( 1,as.integer(dim(df)[1]/mywindow) + 1 ))) {
		currtotalwidth  = totalwidth / $PDFSCALE
		currtotalheight = totalheight_nopk_last / $PDFSCALE
		currtotalratio = totalratio_nopk_last
		pdf(pdfout_nopk,width=currtotalwidth,height=currtotalheight)
		grid.arrange(plot_list[[i]],p2,ncol=1,nrow=mynrow,heights=currtotalratio)
		dev.off()
	} else {
		currtotalwidth = totalwidth / $PDFSCALE
		currtotalheight = totalheight_nopk / $PDFSCALE
		pdf(pdfout_nopk,width=currtotalwidth,height=currtotalheight)
		grid.arrange(plot_list[[i]],ncol=1)
		dev.off()
	}
}

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
=comment
			if (($STRAND eq "Pos" and $strand3 eq "Pos" and $type3 eq "CH") or ($STRAND eq "Neg" and $strand3 eq "Neg" and $type3 eq "GH")) {
				$pngoutDir = $p == 0 ? $PEAK->[0] : $PEAK->[1]; #"PEAK" : "NOPK";
			}
			elsif (($STRAND eq "Pos" and $strand3 eq "Neg" and $type3 eq "GH") or ($STRAND eq "Neg" and $strand3 eq "Pos" and $type3 eq "CH")) {
				$pngoutDir = $p == 0 ? $PEAK->[0] . $TEMP->[1] : $PEAK->[1] . $TEMP->[1]; #"PEAKNEG/" : "NOPKNEG/";
			}
			elsif (($STRAND eq "Pos" and $strand3 eq "Pos" and $type3 eq "CG") or ($STRAND eq "Neg" and $strand3 eq "Neg" and $type3 eq "GC")) {
				$pngoutDir = $p == 0 ? $PEAK->[0] . $CPG->[1] : $PEAK->[1] . $CPG->[1]; #"PEAK/" : "NOPK/"; CG
			}
			elsif (($STRAND eq "Pos" and $strand3 eq "Neg" and $type3 eq "GC") or ($STRAND eq "Neg" and $strand3 eq "Pos" and $type3 eq "CG")) {
				$pngoutDir = $p == 0 ? $PEAK->[0] . $TEMP->[1] . $CPG->[1] : $PEAK->[1] . $TEMP->[1] . $CPG->[1]; #"PEAKNEG/" : "NOPKNEG/"; CG
			}
			else {
				$pngoutDir = "ALL/";
			}
=cut

__END__
	$R->{PDF} = "

# PDF
pdf(\"$pdfout\",width=totalwidth/100,height=totalheight/100)
if (mynrow == 3) {
	grid.arrange(p.pdf,p2.pdf,p3.pdf,ncol=1,nrow=mynrow,heights=totalratio)
} else {
	grid.arrange(p.pdf,p2.pdf,ncol=1,nrow=mynrow,heights=totalratio)
}
dev.off()

# PDF all Conv
pdfout_peak_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_peak_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2.pdf)
dev.off()

";
	


	$R->{PDF_nopk} = "
# PDF
totalheight = (dim(df.rand.1000)[1] + 31.25) * myscale
pdf(\"$pdfout\",width=totalwidth,height=totalheight)
grid.arrange(p.rand,p2.rand,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()
	
# PDF all Conv
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=totalwidth,height=31.25*myscale)
grid.arrange(p2)
dev.off()

";



__END__
	$R->{PDF_nopk} = "

# PDF
currheight = (dim(df.rand.1000)[1] + 31.25) * myscale / 1000
currwidth = totalwidth / 100
print(currheight)
pdf(\"$pdfout\",width=currwidth,height=currheight)
grid.arrange(p.rand.1000.pdf,p2.rand.1000.pdf,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / 100
currwidth = totalwidth / 100
print(currheight)
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.ALL.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=currwidth,height=currheight)
grid.arrange(p2.pdf)
dev.off()

";

	$R->{PDF_nopk_rand_subsample} = "
currheight = (dim(df.rand.subsample)[1]) * myscale / 100
currwidth = totalwidth / 100
pdfout_nopk_rand_subsample = \"$pdfoutFolder/ALL/$pdfoutFilename.rand.subsample.pdf\"
# PDF
pdf(pdfout_nopk_rand_subsample,width=currwidth,height=currheight)
grid.arrange(p.rand.subsample.pdf,p2.rand.subsample.pdf,ncol=1,nrow=mynrow,heights=totalratio)
dev.off()

# PDF all Conv
currheight = 31.25 * myscale / 100
currwidth = totalwidth / 100
pdfout_nopk_all_c_conv = \"$pdfoutFolder/ALL/$pdfoutFilename.rand.subsample.c_conv.pdf\"
pdf(pdfout_nopk_all_c_conv,width=currwidth,height=currheight)
grid.arrange(p2.rand.subsample.pdf)
dev.off()

";




__END__
# BAD PNG TO PDF


__END__

				# Cluster and Third Plot Cluster Graph

				# Read Bed File
#				if (-e $bedFile and linecount($bedFile) > 0 and $currFile !~ /\.NOPK\./) {
					#$Rscript .= $R->{peakbedFile};
#				}

				# Main Plot
#				if ($currFile =~ /\.NOPK\./) {
#					$Rscript .= $R->{mainplot_nopk};
#					$Rscript .= $R->{mainplot_nopk_rand_1000};
#				}
#				else {
#					print "\n\n----------------- $LCY$currFile$N IS A PEAK FILE R = $currFile.PNG.$subsample.R -------------- \n\n";
#					$Rscript .= $R->{mainplot}; # p png and p pdf
#					$Rscript .= $R->{mainplot_nopk_rand_100};
#				}

				# Main Plot Cluster Addition
#				if (-e $curr_cluster_file and -s $curr_cluster_file > 0 and $currFile !~ /\.NOPK\./) {
#					$Rscript .= $R->{mainplotClusterAddition};
#				}
				# Main Plot Peak Bed Addition
#				if (-e $bedFile and linecount($bedFile) > 0 and $currFile !~ /\.NOPK\./) {
#					$Rscript .= $R->{mainplotPeakBedAddition};
#				}
				
				# Second Plot & Conversion Graph
#				if ($currFile !~ /\.NOPK\./) {
#					$Rscript .= $R->{secondplotConversionGraph};
#				}
#				else {
#					$Rscript .= $R->{secondplotConversionGraph};
#					$Rscript .= $R->{secondplotConversionGraph_rand_1000};
#					$Rscript .= $R->{secondplotConversionGraph_rand_100};
#				}

#				if (defined $opt_B and -e $opt_B) {
#					if ($currFile !~ /\.NOPK\./) {
#						LOG($outLog, date() . "\t\t-> ADDED $boxFile!\n","NA") if defined $boxFile;
#						$Rscript .= $R->{box};
#					}
#					else {
#						LOG($outLog, date() . "\t\t-> ADDED $boxFile!\n","NA") if defined $boxFile;
#						$Rscript .= $R->{box_nopk};
#					}
#				}
				# Main Plot
				# Add Third Plot and Do PNG
#				$Rscript .= $R->{Scale};
#				# Main Plot
#				if ($currFile !~ /\.NOPK\./) {
#					$RscriptPNG = $Rscript . $R->{PNG};
#					$RscriptPDF = $Rscript . $R->{PDF};
#				}
#				else {
#					$RscriptPNG = $Rscript . $R->{PNG_nopk};
#					$RscriptPDF = $Rscript . $R->{PDF_nopk};
#					$RscriptPNG .= $Rscript . $R->{PNG_nopk_rand_100};
#					$RscriptPDF .= $Rscript . $R->{PDF_nopk_rand_100};
#					$RscriptPNG_nopk_ALL = $Rscript . $R->{PNG_nopk_ALL};
#					$RscriptPDF_nopk_ALL = $Rscript . $R->{PDF_nopk_ALL};
#				}


