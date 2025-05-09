#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_w $opt_g $opt_G $opt_v $opt_n $opt_r $opt_R $opt_B $opt_c $opt_F $opt_0 $opt_J $opt_i $opt_I $opt_C $opt_0 $opt_f $opt_p);
getopts("n:vg:w:G:r:R:B:cFf0J:i:I:C0p");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/softwares/';
}

use myFootLib;
use FAlite;
use footLoop_sbatch;

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

-0 dry run/debug

-G $LGN<gene to process>$N] 

default: Run ONLY CH/GH
-c: Run ONLY CG/GC (use this for in vitro data)
-C: Run BOTH CG/GC and CH/GH

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
my ($file) = $opt_i;
main($opt_n, $file);

sub main {
	my @Rscript;
	my ($resDir, $file) = @_;
	my ($folder1, $fileName1) = getFilename($file,"folderfull");
	
	my $resDir2 = getFullpath($resDir);
	my ($resDirFullpath) = getFullpath($resDir);
	my %scp;
   makedir("$resDir/.CALL") if not -d "$resDir/\.CALL";
   makedir("$resDir/PEAKS_GENOME") if not -d "$resDir/PEAKS_GENOME";
   makedir("$resDir/PEAKS_LOCAL") if not -d "$resDir/PEAKS_LOCAL";
   makedir("$resDir/PNG") if not -d "$resDir/PNG";
	
	system("mkdir -p $resDir/.footPeak_graph_sbatch_2/") if not -d "$resDir/.footPeak_graph_sbatch_2/";
	my ($OUTDIRS, $PEAK, $TEMP, $RCONV, $CPG, $ALL);
	($OUTDIRS->{PNG}, $PEAK, $TEMP, $RCONV, $CPG, $ALL) = makeOutDir($resDirFullpath . "/PNG/");
	($OUTDIRS->{PDF}) = makeOutDir($resDirFullpath . "/PDF/");
	my %files;
	my ($footPeak_logFile) = "$resDir/footPeak_logFile.txt";
	my $footPeak_graph_logFile = "$resDir/.footPeak_graph_sbatch_2/.$fileName1\_footPeak_graph_sbatch_2_logFile.txt";
	open (my $outLog, ">", $footPeak_graph_logFile) or die "\n\nFailed to write to $footPeak_graph_logFile: $!\n\n";

	LOG($outLog, ">footPeak_graph_sbatch_2.pl version $version\n");
	LOG($outLog, ">UUID: $uuid\n", "NA");
	LOG($outLog, ">Date: $date\n", "NA");
	LOG($outLog, ">Run script: $0 -n $opt_n -i $opt_i\n", "NA");
	
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
- ggplot2 v4.2.0
");

	foreach my $line (@lines) {
		chomp($line);
		if ($line =~ /^[ \t]*def=.+, coor=.+/) {
			$line =~ s/(^\s+|\s+$)//g;
			my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
			LOG($outLog, "gene=$gene,chr=$CHR,beg=$BEG,end=$END,gene=$GENE,val=$VAL,strand=$STRAND\n","NA");
	   	if (defined $opt_G and $gene !~ /$opt_G/i) {
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
      my @strands = qw(Pos Neg);
      for (my $h1 = 0; $h1 < @strands; $h1++) {
         for (my $h2 = 0; $h2 < 4; $h2++) {
            my $strand = $strands[$h1];
            my $type = $types[$h2];
            my $thres2 = $thres; $thres =~ s/0$//;
            my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK";
            $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres2\_$type.PEAK" if not -e $peakFile;
		   	if (defined $opt_G and $peakFile !~ /$opt_G/i) {
		   	   LOG($outLog, date() . " Skipped $LCY$peakFile$N as it doesn't contain $LGN-G $opt_G$N\n");
		   	   next;
		   	}
            LOG($outLog, "label=$label, gene$mygene, strand=$strand, peak=$peakFile\n","NA");
            $files{$peakFile} = $mygene;
         }
      }
   }


	my ($lastGENE, $totalFile, $currfileCount, $lastfile, $fileCount) = (-1,scalar(keys %files),0,-1,0);
	my %Rscripts;

	my $GENE = $files{$file};
	my $STRAND = $coor{$GENE}{STRAND};
	my ($RDSTRAND, $RDCONVTYPE) = $file =~ /^.+_(Pos|Neg|Unk)_.+_(CG|CH|GC|GH).PEAK$/;
	my $RDFLAG = getFlag($file, $STRAND, $RDSTRAND, $RDCONVTYPE);

	my $RDMADEPNG = $opt_r;
	$RDMADEPNG = 0 if $RDFLAG =~ /(RCONV)/ and $opt_r < 2;
	$RDMADEPNG = 0 if $RDFLAG =~ /(ALL)/ and $opt_r < 3;
	$RDMADEPNG = 0 if $RDFLAG =~ /_C$/ and not defined $opt_c and not defined $opt_C;
	$RDMADEPNG = 0 if $RDFLAG !~ /_C$/ and defined $opt_c and not defined $opt_C;

	my $RDMADEPDF = $opt_R;
	$RDMADEPDF = 0 if $RDFLAG =~ /(RCONV)/ and $opt_R < 2;
	$RDMADEPDF = 0 if $RDFLAG =~ /(ALL)/ and $opt_R < 3;
	$RDMADEPDF = 0 if $RDFLAG =~ /_C$/ and not defined $opt_c and not defined $opt_C;
	$RDMADEPDF = 0 if $RDFLAG !~ /_C$/ and defined $opt_c and not defined $opt_C;

	LOG($outLog, "$fileCount $LCY$file$N $LGN$STRAND$N $LCY$RDSTRAND$N $LGN$RDCONVTYPE$N $LPR$RDFLAG$N\n");
	if ($RDMADEPNG eq 0 and $RDMADEPDF eq 0) {
		LOG($outLog, "DONT MAKE MPNG\n");
		return 0;
	}
	else {
		LOG($outLog, "MAKE PNG\n") if $RDMADEPNG ne 0;
		LOG($outLog, "MAKE PDF\n") if $RDMADEPDF ne 0;
	}

	
	my $NA = "NA";
	undef $NA if ($fileCount < 10 or ($fileCount > 10 and $fileCount % 100 == 0));
	LOG($outLog, "\n\n--------------------------\nFILE COUNT = $fileCount\n",$NA);
	my $geneStrandPrint = $STRAND eq "Pos" ? "$LRD$STRAND$N" : $STRAND eq "Neg" ? "$LCY$STRAND$N" : "$LGN$STRAND$N";

  	if (defined $opt_G and $file !~ /$opt_G/i) {
      LOG($outLog, date() . " Skipped $LCY$file$N as it doesn't contain $LGN-G $opt_G$N\n",$NA);
      return 0;
   }
	if ($GENE ne $lastGENE) {
		LOG($outLog, "\n$YW -------- $fileCount/$totalFile Doing$LPR $GENE$N (STRAND = $geneStrandPrint) ---------$N\n\n",$NA) if $GENE ne $lastGENE;
		$currfileCount = 0;
	}
	$lastGENE = $GENE;
	if (defined $opt_g and $GENE ne $opt_g) {
		LOG($outLog, date() . " $LCY Skipped $GENE$N (requested gene is $LGN$opt_g$N\n",$NA);
		return 0;
	} 
	if (defined $opt_w and $file !~ /$opt_w/) {
		LOG($outLog, date() . " Skipped $LCY$file$N (requested want is $LGN$opt_w$N\n",$NA);
		return 0;
	} 

	$lastfile = $file;
	return 0 if not defined $files{$file};

	my ($pk_filename) =  getFilename($file, 'full') . ".out";
	my $peakFile      =  "$resDir/.CALL/$pk_filename";
	my $nopkFile      =  $peakFile;
	my $cluster_file  = "$resDir/FOOTCLUST/CLUST_LOCAL/" . getFilename($file, "full") . ".local.bed.indiv.clust";
	my $kmer_file     =  "$resDir/FOOTCLUST/.TEMP/$pk_filename";
		$nopkFile      =~ s/\.PEAK.out$/.NOPK.out/;
		$kmer_file     =~ s/.out$/.local.bed.clust.kmer/;
	my $bedFile       =  "$resDir/PEAKS_LOCAL/$pk_filename.local.bed";
		$bedFile       =~ s/.out.local.bed/.local.bed/;
	my $totpeak = (-e $peakFile and -s $peakFile > 0) ? linecount($peakFile) : 0;
	my $totnopk = (-e $nopkFile and -s $nopkFile > 0) ? linecount($nopkFile) : 0;
	my $mem = 4000;
	if ($totpeak > 2000 or $totnopk > 2000) {
		$mem = 32000;
	}
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
		my $flag = getFlag($currFile, $geneStrand, $readStrand, $rconvType);
		my $pngoutDir = $flag;
		my $pdfoutDir = $flag;
		LOG($outLog, "\n" . date() . "$LGN$currfileCount.$N $flag $readStrandPrint $rconvTypePrint $LCY$currFile$N\n","NA");
		LOG($outLog, "$readStrandPrint\t$rconvTypePrint\t$flag\n","NA");



		LOG($outLog, "flag=$LPR $flag$N opt_r=$opt_r\n");

		my $madePNG = $opt_r;
		$madePNG = 0 if $flag =~ /(RCONV)/ and $opt_r < 2;
		$madePNG = 0 if $flag =~ /(ALL)/ and $opt_r < 3;
		$madePNG = 0 if $flag =~ /_C$/ and not defined $opt_c and not defined $opt_C;
		$madePNG = 0 if $flag !~ /_C$/ and defined $opt_c and not defined $opt_C;

		my $madePDF = $opt_R;
		$madePDF = 0 if $flag =~ /(RCONV)/ and $opt_R < 2;
		$madePDF = 0 if $flag =~ /(ALL)/ and $opt_R < 3;
		$madePDF = 0 if $flag =~ /_C$/ and not defined $opt_c and not defined $opt_C;
		$madePDF = 0 if $flag !~ /_C$/ and defined $opt_c and not defined $opt_C;

		$madePNG = 1 if $madePNG > 0;
		$madePDF = 1 if $madePDF > 0;
		my $currFilename2 = $currFilename;
		$currFilename2 =~ s/^.+_bismark_bt2.bam_gene(.+)$/$1/;
		$currFilename2 =~ s/.(NOPK|PEAK).out//;
		my $pngout = "$resDir/PNG/$pngoutDir/$currFilename2.$flag.png";
		my $pngout_conly = "$resDir/PNG/$pngoutDir/CONLY/$currFilename2.$flag.png.conly.png";
		my $pdfout = "$resDir/PDF/$pdfoutDir/$currFilename2.$flag.pdf";
		my $pdfout_conly = "$resDir/PDF/$pdfoutDir/CONLY/$currFilename2.$flag.pdf.conly.pdf";
		my $lenpdfout = "$resDir/PNG/$pngoutDir/LENGTH/$currFilename2\_length.pdf";
		LOG($outLog, date() . " --> DEBUG flag=$flag madePNG = $madePNG pngout=$pngout\n");

		my ($RscriptPDF, $RscriptPNG, $RscriptPDF_nopk_ALL, $RscriptPNG_nopk_ALL);

		my $Rscript = "
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)
";
		
		my $totread = $totpeak + $totnopk;
		if (not -e $currFile or (-e $currFile and -s $currFile <= 2)) { 
			$Rscript .= "png(\"$pngout\",1000,1000)\nplot(NA,xlim=c(1,100),ylim=c(1,100),xlab=NA,ylab=NA,bty=\"n\")\ntext(50,50,cex=3,labels=c(\"$currFilename2\n\nPEAK = $totpeak / $totread\"))\ndev.off()\n";
			$Rscript .= "png(\"$pngout_conly\",1000,1000)\nplot(NA,xlim=c(1,100),ylim=c(1,100),xlab=NA,ylab=NA,bty=\"n\")\ntext(50,50,cex=3,labels=c(\"$currFilename2\n\nPEAK = $totpeak / $totread\"))\ndev.off()\n";
			$RscriptPNG = $Rscript;
			$RscriptPDF = $Rscript;
			$RscriptPNG_nopk_ALL = $Rscript;
			$RscriptPDF_nopk_ALL = $Rscript;
			LOG($outLog, "${LPR}PNGOUT$N:\n$LCY$pngout$N $LGN #Not run coz totpeak is $totpeak$N\n\n");
			LOG($outLog, "${LPR}PDFOUT$N:\n$LCY$pdfout$N $LGN #Not run coz totpeak is $totpeak$N\n\n");
		}
		else {
			open (my $incurrFile, "cut -f1 $currFile|") or DIELOG($outLog, "Failed to open $currFile: $!\n");
			open (my $outcurrFileID, ">", "$currFile.id") or DIELOG($outLog, "Failed to open $currFile.id: $!\n");
			while (my $line = <$incurrFile>) {
				chomp($line);
				my ($num1, $num2, $num3) = $line =~ /^.*m([A-Za-z0-9]+_\d+)_.+\/(\d+)\/(ccs|\d+_\d+)/;
				DIELOG($outLog, "\n\n" . date() . " Failed to parse numbers from currFile=$LCY$currFile$N id = $LPR$line$N\n\n") if not defined $num1 or not defined $num2 or not defined $num3;
				$num3 = 0 if $num3 eq "ccs";
				my $num = "$num1$num2$num3";
		      $num =~ s/_//g;
				print $outcurrFileID "$line\t$num\n";
			}
			close $incurrFile;
			close $outcurrFileID;
			
			my ($strand, $convtype) = $currFile =~ /(Pos|Neg).+(CG|CH|GC|GH).(PEAK|NOPK).+$/;
			my $parseName = parseName($currFile);
			my $plasmid = $parseName->{plasmid};
			my $desc = $parseName->{desc};
			my $bc = $parseName->{bc};
			my $label = $parseName->{label}; $label =~ s/_bismark_bt2.bam//;
			my $mytitle = "$label BC$bc $plasmid $desc $strand $convtype $flag";
			my ($R) = Rscript($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout, $resDir, $convtype, $mytitle);

			# Read Table
			$Rscript .= $R->{readTable};

			# Cluster and Third Plot Cluster Graph
			if (-e $curr_cluster_file and -s $curr_cluster_file > 0 and $currFile !~ /\.NOPK\./) {
				$Rscript .= $R->{clusterFile};
			} 
			else {
				$Rscript .= $R->{noclusterFile};
			}

			# Read Bed File
			if (-e $bedFile and -s $bedFile > 0 and $currFile !~ /\.NOPK\./) {
				$Rscript .= $R->{peakbedFile};
			}

			# Main Plot
			if ($currFile =~ /\.NOPK\./) {
				$Rscript .= $R->{mainplot_nopk}; 
				$Rscript .= $R->{mainplot_nopk_rand_1000};
				$Rscript .= $R->{mainplot_nopk_rand_100}; 
			}
			else {
				$Rscript .= $R->{mainplot}; 
				$Rscript .= $R->{mainplot_peak_rand_300}; 
			}

			# Main Plot Cluster Addition
			if (-e $curr_cluster_file and -s $curr_cluster_file > 0 and $currFile !~ /\.NOPK\./) {
				$Rscript .= $R->{mainplotClusterAddition};
			}
			# Main Plot Peak Bed Addition
			if (-e $bedFile and -s $bedFile > 0 and $currFile !~ /\.NOPK\./) {
				$Rscript .= $R->{mainplotPeakBedAddition};
			}
			
			# Second Plot & Conversion Graph
			if ($currFile !~ /\.NOPK\./) {
				$Rscript .= $R->{secondplotConversionGraph};
			}
			else {
				$Rscript .= $R->{secondplotConversionGraph};
				$Rscript .= $R->{secondplotConversionGraph_rand_1000};
				$Rscript .= $R->{secondplotConversionGraph_rand_100};
			}

			if (defined $opt_B and -e $opt_B) {

				my ($pk_filename_short) = $pk_filename =~ /^.+gene(.+)_(Pos|Neg|Unk)_\d+.+$/;
				my ($pk_filename_plasmid) = $pk_filename_short =~ /^(.+)_desc.+$/i;
				$pk_filename_plasmid = $pk_filename_short if not defined $pk_filename_plasmid;
				if ($currFile !~ /\.NOPK\./) { #PEAK
					LOG($outLog, date() . "\t\t-> ADDED $boxFile!\n","NA") if defined $boxFile;
					my ($Rbox_peak, $Rbox_nopk) = getbox($boxFile,$pk_filename_short,$pk_filename_plasmid);
					$Rscript .= $Rbox_peak;
				}
				else { #NOPK
					LOG($outLog, date() . "\t\t-> ADDED $boxFile!\n","NA") if defined $boxFile;
					my ($Rbox, $Rbox_nopk) = getbox($boxFile,$pk_filename_short,$pk_filename_plasmid);
					$Rscript .= $Rbox_nopk;
				}
			}

			# Add Third Plot and Do PNG
			$Rscript .= $R->{Scale};
			# Main Plot
			$RscriptPNG = $Rscript;
			$RscriptPDF = $Rscript;
			$RscriptPNG_nopk_ALL = $Rscript;
			$RscriptPDF_nopk_ALL = $Rscript;

			if ($currFile !~ /\.NOPK\./) {
				$RscriptPNG .= $R->{PNG};
				$RscriptPNG .= $R->{PNG_peak_300};
				$RscriptPDF .= $R->{PDF};
				$RscriptPDF .= $R->{PDF_peak_300};
			}
			else {
				$RscriptPNG .= $R->{PNG_nopk};
				$RscriptPNG .= $R->{PNG_nopk_rand_100};
				$RscriptPNG_nopk_ALL .= $R->{PNG_nopk_ALL};

				$RscriptPDF .= $R->{PDF_nopk};
				$RscriptPDF .= $R->{PDF_nopk_rand_100};
				$RscriptPDF_nopk_ALL .= $R->{PDF_nopk_ALL};
			}
		}

		open (my $outRscriptPNG, ">", "$currFile.PNG.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PNG.R: $!\n") and print $outLog $Rscript and next);
		print $outRscriptPNG $RscriptPNG;
		$Rscripts{"$currFile.PNG.R"}{graphfile} = $pngout;
		$Rscripts{"$currFile.PNG.R"}{summary} = $summary;
		$Rscripts{"$currFile.PNG.R"}{runR} = $madePNG; 
		$Rscripts{"$currFile.PNG.R"}{runType} = $flag;
		close $outRscriptPNG;
		LOG($outLog, "${LPR}PNGOUT$N:\n$LCY$pngout$N # $LGN$totpeak$N peak\n\n") if $madePNG ne 0;

		open (my $outRscriptPDF, ">", "$currFile.PDF.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PDF.R: $!\n") and print $outLog $Rscript and next);
		print $outRscriptPDF $RscriptPDF;
		close $outRscriptPDF;
		$Rscripts{"$currFile.PDF.R"}{graphfile} = $pdfout;
		$Rscripts{"$currFile.PDF.R"}{summary} = $summary;
		$Rscripts{"$currFile.PDF.R"}{runR} = $madePDF; 
		$Rscripts{"$currFile.PDF.R"}{runType} = $flag;
		LOG($outLog, "${LPR}PDFOUT$N:\n$LCY$pdfout$N # $LGN$totpeak$N peak\n\n") if $madePDF ne 0;

		if ($currFile =~ /\.NOPK\./) {
			open (my $outRscriptPNG_nopk_ALL, ">", "$currFile.PNG_nopk_ALL.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PNG_nopk_ALL.R: $!\n") and print $outLog $Rscript and next);
			print $outRscriptPNG_nopk_ALL $RscriptPNG_nopk_ALL;
			$Rscripts{"$currFile.PNG_nopk_ALL.R"}{summary} = $summary;
			$Rscripts{"$currFile.PNG_nopk_ALL.R"}{runR} = 0;
			$Rscripts{"$currFile.PNG_nopk_ALL.R"}{runType} = $flag . "_ALL";
			close $outRscriptPNG_nopk_ALL;
		}



		if ($currFile =~ /\.NOPK\./) {
			open (my $outRscriptPDF_nopk_ALL, ">", "$currFile.PDF_nopk_ALL.R") or (LOG($outLog, date() . "Failed to write R script into $currFile.PDF_nopk_ALL.R: $!\n") and print $outLog $Rscript and next);
			print $outRscriptPDF_nopk_ALL $RscriptPDF_nopk_ALL;
			$Rscripts{"$currFile.PDF_nopk_ALL.R"}{summary} = $summary;
			$Rscripts{"$currFile.PDF_nopk_ALL.R"}{runR} = 0;
			$Rscripts{"$currFile.PDF_nopk_ALL.R"}{runType} = $flag . "_ALL";
			close $outRscriptPDF_nopk_ALL;
		}
		LOG($outLog, "${LPR}PDFOUT$N:\n$LCY$pdfout$N # $LGN$totpeak$N peak\n\n");
	}
	$totalFile = (keys %Rscripts);
	LOG($outLog, "\n\n$YW ----------------- Running $totalFile/$fileCount R Scripts (below, showing only that are run) ------------------$N\n\n");
	# open outRscripts for Rscripts that aren't relevant
	$fileCount = 0;
	foreach my $outRscriptGRAPH (sort keys %Rscripts) {
		my $runR = $Rscripts{$outRscriptGRAPH}{runR};
		my $summary = $Rscripts{$outRscriptGRAPH}{summary};
		my $pngout = $Rscripts{$outRscriptGRAPH}{graphfile};
		LOG($outLog, "${LPR}Rscript$N: $YW$outRscriptGRAPH$N\n${LPR}GRAPH    $N: $LCY$pngout$N\n${LPR}Summary$N: $LGN$summary$N\n\n") if $runR == 1;
		push(@Rscript, $outRscriptGRAPH) if $runR == 1;
	}
	
	$max_parallel_run = scalar(@Rscript);
	my $sbatch_cmd = "Rscript FILENAME";
	my $force_sbatch = 1 if defined $opt_F;
	my $outsbatchDir = "$resDir/.footPeak_graph_sbatch_2/";
	system("mkdir -p $outsbatchDir") if not -d $outsbatchDir;
	my $debug;

	my $forcerun = "off";
	   $forcerun = "on" if defined $opt_F;
	   $forcerun = "on" if defined $opt_f;
	LOG($outLog, "forcerun is $YW$forcerun$N\n\n");

	if (defined $opt_p) {
	   LOG($outLog, "#Parallel run with slurm enabled (-p)\n\n");
	   if (defined $opt_0) {
	      LOG($outLog, date() . "\n" . "$sbatch_cmd\n");
	   }
	   my $mem = 16000;
	   my $debug = 1 if defined $opt_0;
	   my ($max_parallel_run) = defined $opt_J ? $opt_J : 1;
	   my $force_sbatch = $opt_F;
  	   $force_sbatch = $opt_f if not defined $opt_F;
		footLoop_sbatch_main($sbatch_cmd, "footPeak_graph_sbatch_2", \@Rscript, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir, $mem, $debug);
	}
	else {
   	LOG($outLog, "#Single run\n\n");
   	for (my $i = 0; $i < @Rscript; $i++) {
	      my $Rscriptindiv = $Rscript[$i];
	      my ($Rscriptindivname) = getFilename($Rscriptindiv, "full");
	      my $RscriptindivDone = "$outsbatchDir/$Rscriptindivname.done";
	      my $currcmd = $sbatch_cmd;
	      my $indice = $i + 1;
	      $currcmd =~ s/FILENAME/$Rscriptindiv/g;
		      $currcmd =~ s/FNINDICE/$indice/g;
	      if (defined $opt_0) {
  		      LOG($outLog, date() . "$YW$indice/$totalFile$N $LCY$Rscriptindivname$N\n");
   	      LOG($outLog, "\t$LGN$currcmd$N\n");
   	   }
   	   elsif (not -e $RscriptindivDone or defined $opt_F or defined $opt_f) {
				LOG($outLog, "\n\t$LCY$currcmd$N\n");
   	      system("$currcmd > $outsbatchDir/$Rscriptindivname.log.txt 2>&1") == 0 or LOG($outLog, "Failed to run:\n$LCY$currcmd$N\n$outsbatchDir/$Rscriptindivname.log.txt 2>&1\n$LRD$!$N\n\n");
   	      system("touch $RscriptindivDone") if not -e $RscriptindivDone;
   	   }
   	   else {
   	      LOG($outLog, date() . "${LPR}footPeak.pl$N: $YW$indice/$totalFile$N $LCY$Rscriptindivname$N: using previously made peaks\n");
   	      LOG($outLog, "Done=$LCY$RscriptindivDone$N)\n","NA");
			}
      }
   }
}

###############
# Subroutines #
###############

sub Rscript {
	my ($currFile, $bedFile, $curr_cluster_file, $totpeak, $totnopk, $pngout, $pdfout, $lenpdfout, $resDir, $convtype, $mytitle) = @_;
	my $currFileID = $currFile . ".id";
	my $R;
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
	my $currFilename2 = $currFilename;
	$currFilename2 =~ s/^(.+)_bismark_bt2.bam_/$1_/;
	$currFilename2 =~ s/.(NOPK|PEAK).out//;

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
p2.png.scale = p.png.scale
p.pdf.scale = 0.33
p2.pdf.scale = p.pdf.scale

#####################
# Read Table
convwant = \"$convtype\"
peakFile = gsub(\".out\$\",\"\",\"$currFile\")

myseq = read.table(peakFile,sep=\"\\t\",nrow=1,colClasses = c(\"character\"))
myseq = myseq[,seq(6,dim(myseq)[2])]
myseq = data.frame(variable=seq(1,dim(myseq)[2]),value=t(myseq))
colnames(myseq) = c(\"variable\",\"value\")
myseq\$value0 = c(\"N\",myseq\$value[1:(dim(myseq)[1]-1)])
myseq\$value2 = c(myseq\$value[2:(dim(myseq)[1])],\"N\")

CGprof = function(df,window=200,step=10) {
	df\$value0 = c(\"N\",df\$value[1:(dim(df)[1]-1)])
	df\$value2 = c(df\$value[2:(dim(df)[1])],\"N\")
	CpG = c()
	GC = c()
	GCskew = c()
	ATskew = c()
	ind = 0
	for (i in seq(1,dim(df)[1],10)) {
		ind = ind + 1
		i0 = max(i - window/2,1)
		i1 = min(i + window/2,dim(df)[1])
		temp = df[df\$variable >= i0 & df\$variable <= i1,]
		nC = dim(temp[temp\$value == \"C\",])[1]
		nG = dim(temp[temp\$value == \"G\",])[1]
		nT = dim(temp[temp\$value == \"T\",])[1]
		nA = dim(temp[temp\$value == \"A\",])[1]
		nCG = dim(temp[temp\$value == \"C\" & temp\$value2 == \"G\",])[1]
		nTot = nA+nG+nT+nC
		if (nTot > 0) {
			CpG[ind] = (nCG * nTot) / (nG * nC) 
			GC[ind] = (nC + nG)/(nA+nG+nT+nC)
			if (nG+nC > 0) {
				GCskew[ind] = (nG-nC)/(nG+nC)
			} else {
				GCskew[ind] = 0
			}
			if (nA+nT > 0) {
				ATskew[ind] = (nA-nT)/(nA+nT)
			} else {
				ATskew[ind] = 0
			}
		} else {
			CpG[ind] = 0
			GC[ind] = 0
			GCskew[ind] = 0
			ATskew[ind] = 0
		}
		
	}
	mylist = list(x=seq(1,dim(df)[1],10),CGdens=CpG,GCcont=GC,GCskew=GCskew,ATskew=ATskew)
	dm = data.frame(x=seq(1,dim(df)[1],10),CGdens=CpG,GCcont=GC,GCskew=GCskew,ATskew=ATskew)
	return(invisible(dm))
}

myres = CGprof(subset(myseq,select=c(\"variable\",\"value\")))
myres\$CGdens = (myres\$CGdens - 0) / 2 * 1
myres\$GCcont = (myres\$GCcont - 0) / 1 * 1
myres\$GCskew = (myres\$GCskew + 1) / 2 * 1
myres\$ATskew = (myres\$ATskew + 1) / 2 * 1
library(ggplot2)
library(reshape2)
myresdm = melt(myres,id.vars=\"x\")
p3.CGprof.col = c(\"CGdens\"=\"blue2\",\"GCcont\"=\"green4\",\"GCskew\"=\"red3\",\"ATskew\"=\"orange\")
p3.CGprof = ggplot(myresdm,aes(x,value)) +
	geom_line(aes(color=variable),lwd=0.1) +
	ylim(c(0,1)) +
	theme_bw() + theme(panel.grid = element_blank()) +
	ylab(\"\") + xlab(\"bp\") +
	scale_color_manual(values=p3.CGprof.col,expand=c(0,0))

p3.CGprof.mod.png = ggplot(myresdm,aes(x,value)) +
	geom_line(aes(color=variable),lwd=0.1) +
	ylim(c(0,1)) +
	theme_bw() + theme(panel.grid = element_blank()) +
	ylab(\"\") + xlab(\"bp\") +
	scale_color_manual(values=p3.CGprof.col,expand=c(0,0))


p3.CGprof.mod.pdf = ggplot(myresdm,aes(x,value)) +
	geom_line(aes(color=variable),lwd=0.1) +
	ylim(c(0,1)) +
	theme_bw() + theme(panel.grid = element_blank()) +
	ylab(\"\") + xlab(\"bp\") +
	scale_color_manual(values=p3.CGprof.col,expand=c(0,0))

p3.CGprof.mod.png = p3.CGprof.mod.png +
	annotate(geom=\"segment\",x=0.04*max(myresdm\$x),xend=0.04*max(myresdm\$x),y=0.0,yend=1.00,color=\"blue2\") +
	annotate(geom=\"segment\",x=0.08*max(myresdm\$x),xend=0.08*max(myresdm\$x),y=0.0,yend=1.00,color=\"red3\")
for (i in seq(0,1,0.25)) {
	p3.CGprof.mod.png = p3.CGprof.mod.png + annotate(geom=\"segment\",x=0.03*max(myresdm\$x),xend=0.04*max(myresdm\$x),y=i,yend=i,color=\"blue2\")
	p3.CGprof.mod.png = p3.CGprof.mod.png + annotate(geom=\"text\"   ,x=0.03*max(myresdm\$x),y=i,label=i*2,size=5*p2.png.scale,hjust=1,color=\"blue2\")
	p3.CGprof.mod.png = p3.CGprof.mod.png + annotate(geom=\"segment\",x=0.08*max(myresdm\$x),xend=0.09*max(myresdm\$x),y=i,yend=i,color=\"red3\")
	p3.CGprof.mod.png = p3.CGprof.mod.png + annotate(geom=\"text\"   ,x=0.09*max(myresdm\$x),y=i,label=(i*2)-1,size=5*p2.png.scale,hjust=0,color=\"red3\")
	p3.CGprof.mod.png = p3.CGprof.mod.png + annotate(geom=\"text\"   ,x=0.06*max(myresdm\$x),y=i,label=i,size=5*p2.png.scale,hjust=0.5,color=\"green4\")
}
p3.CGprof.mod.png = p3.CGprof.mod.png +
	scale_x_continuous(expand = c(0,0)) +
	theme_blank

p3.CGprof.mod.pdf = p3.CGprof.mod.pdf +
	annotate(geom=\"segment\",x=0.04*max(myresdm\$x),xend=0.04*max(myresdm\$x),y=0.0,yend=1.00,color=\"blue2\",lwd=0.1) +
	annotate(geom=\"segment\",x=0.08*max(myresdm\$x),xend=0.08*max(myresdm\$x),y=0.0,yend=1.00,color=\"red3\",lwd=0.1)
for (i in seq(0,1,0.25)) {
	p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + annotate(geom=\"segment\",x=0.03*max(myresdm\$x),xend=0.04*max(myresdm\$x),y=i,yend=i,color=\"blue2\",lwd=0.1)
	p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + annotate(geom=\"text\"   ,x=0.03*max(myresdm\$x),y=i,label=i*2,size=1,hjust=1,color=\"blue2\")
	p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + annotate(geom=\"segment\",x=0.08*max(myresdm\$x),xend=0.09*max(myresdm\$x),y=i,yend=i,color=\"red3\",lwd=0.1)
	p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + annotate(geom=\"text\"   ,x=0.09*max(myresdm\$x),y=i,label=(i*2)-1,size=1,hjust=0,color=\"red3\")
	p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + annotate(geom=\"text\"   ,x=0.06*max(myresdm\$x),y=i,label=i,size=1,hjust=0.5,color=\"green4\")
}
p3.CGprof.mod.pdf = p3.CGprof.mod.pdf +
	scale_x_continuous(expand = c(0,0)) +
	theme_blank

myseq\$type = \"N\"
myseq[myseq\$value == \"C\",]\$type = \"CH\"
myseq[myseq\$value == \"G\",]\$type = \"GH\"
myseq[myseq\$value == \"C\" & myseq\$value2 == \"G\",]\$type = \"CG\"
myseq[myseq\$value == \"G\" & myseq\$value0 == \"C\",]\$type = \"GC\"


df = read.table(\"$currFile\",sep=\"\\t\")
if (length(grep(\"^[0-9]+\$\",df[,2])) == 0) {
   df = df[,-1]
}
df.id = read.table(\"$currFileID\",sep=\"\\t\",colClasses=c(\"factor\",\"factor\"))
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
		ymax=ymax,x=xmin,y=ymin,color=as.factor(clust2)),
		size=1*p3.png.scale,fill=rgb(0,0,0,alpha=0)) + # SIZE
	geom_rect(data=df4,aes(fill=as.factor(clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax),size=0.5*p3.png.scale) + # SIZE
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
		ymax=ymax,x=xmin,y=ymin,color=as.factor(clust2)),
		size=1*p3.pdf.scale,fill=rgb(0,0,0,alpha=0)) + # SIZE
	geom_rect(data=df4,aes(fill=as.factor(clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax),size=0.5*p3.pdf.scale) + # SIZE
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
if (convwant == \"CH\") {
	dmc = myseq[myseq\$value == \"C\" & myseq\$value2 != \"G\",]
} else if (convwant == \"CG\") {
	dmc = myseq[myseq\$value == \"C\",]
} else if (convwant == \"GH\") {
	dmc = myseq[myseq\$value == \"G\" & myseq\$value0 != \"C\",]
} else if (convwant == \"GC\") {
	dmc = myseq[myseq\$value == \"G\",]
}
dmc = dmc[order(dmc\$variable),]
dmc\$x= seq(1,dim(dmc)[1])
dmc = dmc[order(dmc\$variable),]
dmc = subset(dmc,select=c(variable,x))

i0 = dmc[dim(dmc)[1],,drop=F]
dmc2 = data.frame()
for (i in seq(max(dmc\$variable),1,-1)) {
   if (dim(dmc[dmc\$variable == i,])[1] == 0) {
      icurr = i0
      icurr\$variable = i
      dmc2 = rbind(dmc2,icurr)
   } else {
      i0 = dmc[dmc\$variable == i,,drop=F]
      dmc2 = rbind(dmc2,i0)
   }
}
dmc2 = dmc2[order(dmc2\$variable),]
rownames(dmc2) = seq(1,dim(dmc2)[1])
dmc2\$beg = dmc2\$x
dmc2\$end = dmc2\$x
dmc2\$V2 = dmc2\$variable
dmc2\$V3 = dmc2\$variable


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
totalpeak = $totpeak + $totnopk
if (totalpeak == 0) {
	peakperc = 0
} else {
	peakperc = as.integer(1000*$totpeak / ($totpeak + $totnopk)+0.5)/10
}
p = ggplot(dm,aes(variable,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 theme_bw() + theme(legend.position=\"none\") + 
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL) + theme_blank +
	 ggtitle(paste(\"$mytitle\",\"\\n\",peakperc,\"\% $totpeak/\",totalpeak,sep=\"\"))


p.heatmaponly = ggplot(dm,aes(variable,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL) +
	 theme_void() +
	 theme(legend.position=\"none\")


dm.conly = dm[dm\$variable \%in\% dmc\$variable,]
dm.conly = merge(dm.conly,dmc,by=\"variable\")
p.conly = ggplot(dm.conly,aes(x,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL) +
	 theme_void() +
	 theme(legend.position=\"none\")
p.png = p
p.pdf = p
p.conly.png = p.conly

";


# -------------------- $R->{mainplot_peak_rand_300} 240716
$R->{mainplot_peak_rand_300} .= "
	
	##########################
	# Main Plot Part 1 peak rand
	count0 = apply(df[,-1], 1, function(x){return(length(x[x == 0])/length(x))})
   if (length(count0[count0 >= 0.2]) > 0) {
      print(\"###########################\")
      print(paste(\"Has indel more than 20% of read length:\",length(count0[count0 > 0.2])))
      print(paste(\"Example:\"))
      print(head(df[count0 >= 0.2,1:5]))
      print(\"###########################\")
   }
   if (length(count0[count0 < 0.2]) > 0) {
		df.rand.300 = df[count0 < 0.2,]
	} else {
		print(\"Using all df.rand.300 coz there's almost no reads left\")
		df.rand.300 = df
	}

	df.total = dim(df.rand.300)[1]
	set.seed(42)
	if (df.total > 300) {
		df.rand.300 = df.rand.300[sample( seq(1,df.total) ,300,replace=F),]
		print(\"randoming 300. DF TOTAL = \"); print(df.total)
	} else {
		df.rand.300 = df.rand.300
		print(\"NOT randoming 300. DF TOTAL = \"); print(df.total)
	}

	df.rand.300 = df.rand.300[order(df.rand.300\$y),]
	df.rand.300\$y = seq(1,dim(df.rand.300)[1])

	write.table(df.rand.300,file=\"$currFile.rand\",quote=F,row.names=F,col.names=F,sep=\"\\t\")
	print(\"Wrote to $currFile.rand\")
	
	dm.rand.300 = melt(df.rand.300,id.vars=c(\"V1\",\"y\"))
	print(dim(df.rand.300))
	dm.rand.300\$variable = as.numeric(as.character(dm.rand.300\$variable))
	
	if (convwant == \"CH\") {
		dmc.rand.300 = myseq[myseq\$value == \"C\" & myseq\$value2 != \"G\",]
	} else if (convwant == \"CG\") {
		dmc.rand.300 = myseq[myseq\$value == \"C\",]
	} else if (convwant == \"GH\") {
		dmc.rand.300 = myseq[myseq\$value == \"G\" & myseq\$value0 != \"C\",]
	} else if (convwant == \"GC\") {
		dmc.rand.300 = myseq[myseq\$value == \"G\",]
	}
	dmc.rand.300 = dmc.rand.300[order(dmc.rand.300\$variable),]
	dmc.rand.300\$x= seq(1,dim(dmc.rand.300)[1])
	dmc.rand.300 = subset(dmc.rand.300,select=c(variable,x))
	dmc.rand.300 = dmc.rand.300[order(dmc.rand.300\$variable),]
	i0 = dmc.rand.300[dim(dmc.rand.300)[1],,drop=F]
	dmc2.rand.300 = data.frame()
	for (i in seq(max(dmc.rand.300\$variable),1,-1)) {
	   if (dim(dmc.rand.300[dmc.rand.300\$variable == i,])[1] == 0) {
	      icurr = i0
	      icurr\$variable = i
	      dmc2.rand.300 = rbind(dmc2.rand.300,icurr)
	   } else {
	      i0 = dmc.rand.300[dmc.rand.300\$variable == i,,drop=F]
	      dmc2.rand.300 = rbind(dmc2.rand.300,i0)
	   }
	}
	dmc2.rand.300 = dmc2.rand.300[order(dmc2.rand.300\$variable),]
	rownames(dmc2.rand.300) = seq(1,dim(dmc2.rand.300)[1])
	dmc2.rand.300\$beg = dmc2.rand.300\$x
	dmc2.rand.300\$end = dmc2.rand.300\$x
	dmc2.rand.300\$V2 = dmc2.rand.300\$variable
	dmc2.rand.300\$V3 = dmc2.rand.300\$variable
	
	totalpeak = $totpeak + $totnopk
	if (totalpeak == 0) {
		peakperc = 0
	} else {
		peakperc = as.integer(1000*$totpeak / ($totpeak + $totnopk)+0.5)/10
	}
	p.rand.300 = ggplot(dm.rand.300,aes(variable,y)) +  
		 geom_tile(aes(fill=as.factor(value))) + 
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0)) +
		 theme_void() +
		 theme(legend.position=\"none\") + 
		 ggtitle(paste(\"$mytitle\",\"\\n\",peakperc,\"\% $totpeak/\",totalpeak,sep=\"\"))
	
	p.heatmaponly.rand.300 = ggplot(dm.rand.300,aes(variable,y)) +  
		 geom_tile(aes(fill=as.factor(value))) + 
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0)) +
		 theme_void() +
		 theme(legend.position=\"none\")
	
	dm.rand.300.conly = dm.rand.300[dm.rand.300\$variable \%in\% dmc.rand.300\$variable,]
	dm.rand.300.conly = merge(dm.rand.300.conly,dmc.rand.300,by=\"variable\")
	
	p.rand.300.conly = ggplot(dm.rand.300.conly,aes(x,y)) +  
		geom_tile(aes(fill=as.factor(value))) +
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0)) +
		 labs(x=NULL,y=NULL) +
		 theme_void() +
		 theme(legend.position=\"none\")

	p.rand.300.pdf = p.rand.300 + theme(plot.title = element_text(size = 20*p.pdf.scale))
	p.rand.300.png = p.rand.300 + theme(plot.title = element_text(size = 20*p.png.scale))
	p.rand.300.conly.png = p.rand.300.conly

";


sub getbox {
	my ($boxFile, $boxchr, $boxplasmid) = @_;
	my $Rbox = "
	box = read.table(\"$boxFile\",sep=\"\\t\")
	if ((dim(box[grep(\"^$boxchr\$\",box\$V1,ignore.case=TRUE),])[1] > 0) | (dim(box[grep(\"^$boxplasmid\$\",box\$V1,ignore.case=TRUE),])[1] > 0)) {
		p.col2 = p.col
		if ((dim(box[grep(\"^$boxchr\$\",box\$V1,ignore.case=TRUE),])[1] > 0)) {
			box = box[grep(\"^$boxchr\$\",box\$V1,ignore.case=TRUE),]
         print(paste(\"1 \",dim(box)[1],sep=\"\"))
		} else {
			box = box[grep(\"^$boxplasmid\$\",box\$V1,ignore.case=TRUE),]
         print(paste(\"2 \",dim(box)[1],sep=\"\"))
		}
		if (dim(box)[2] == 7) {
			boxcolordf = as.data.frame(aggregate(box\$V2,by=list(box\$V4,box\$V7),sum))
			colnames(boxcolordf) = c(\"V4\",\"V7\",\"sumz\")
			boxcolor = boxcolordf\$V7
			names(boxcolor) = boxcolordf\$V4
			p.col2 = c(p.col,boxcolor)
			names(p.col2) = c(names(p.col),names(boxcolor))
		} else {
			if (dim(box)[2] <= 3) {
				box\$V4 = box\$V2
			}
			if (dim(box)[2] < 6) {
				box\$V6 = \"+\"
			}
			if (dim(box)[2] < 7) {
				boxcolorname = unique(box\$V4)
				boxcolor = rep(\"#969696\",10)
				while (length(boxcolor) < length(boxcolorname)) {
					boxcolor = c(boxcolor,boxcolor)
				}
				boxcolor = boxcolor[seq(1,length(boxcolorname))]
				names(boxcolor) = boxcolorname
				p.col2 = c(p.col,boxcolor)
				names(p.col2) = c(names(p.col),names(boxcolor))
			}
		}
		box\$label = paste(box\$V4, \" (\",box\$V6,\")\",sep=\"\")
      p.col2[p.col2 == \"#FFFF33\"] = \"black\"
		if (length(p.col2[grep(\"T7.*Prom\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"T7.*Prom\",names(p.col2),ignore.case=TRUE)] = \"#e7298a\"
		}
		if (length(p.col2[grep(\"T3.*Prom\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"T3.*Prom\",names(p.col2),ignore.case=TRUE)] = \"#e7298a\"
		}
		if (length(p.col2[grep(\"Prom\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"Prom\",names(p.col2),ignore.case=TRUE)] = \"#e7298a\"
		}
		if (length(p.col2[grep(\"VR_?[0-9]\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"VR_?[0-9]\",names(p.col2),ignore.case=TRUE)] = \"#045a8d\"
		}
		if (length(p.col2[grep(\"(SNRPN|AIRN)\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"(SNRPN|AIRN)\",names(p.col2),ignore.case=TRUE)] = \"black\"
		}
		if (length(p.col2[grep(\"Barcode\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"Barcode\",names(p.col2),ignore.case=TRUE)] = \"grey\"
		}
		if (length(p.col2[grep(\"Primer\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"Primer\",names(p.col2),ignore.case=TRUE)] = \"grey\"
		}
		if (length(p.col2[grep(\"(sg[0-9]+|sgRNA|Nick)\",names(p.col2),ignore.case=TRUE)]) > 0) {
			p.col2[grep(\"(sg[0-9]+|sgRNA|Nick)\",names(p.col2),ignore.case=TRUE)] = \"orange4\"
		}
	";
	my $Rbox_nopk = $Rbox;


	my $Rbox_peak = $Rbox;
	$Rbox_peak .= "
		print(\"!!! PEAK\")
		print(names(p.col2))
		print(p.col2)
		box2 = box
		box2 = merge(box2,subset(dmc2,select=c(V2,beg)),by=c(\"V2\"))
		box2 = merge(box2,subset(dmc2,select=c(V3,end)),by=c(\"V3\"))
      box2\$x = 0
      box2\$y = 0
      box2\$variable = 0
      box\$x = 0
      box\$y = 0
      box\$variable = 0
		box\$x2 = 0
		box\$y2 = 0
		box\$value = 0

		# p
		p.png = p.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,size=0.5*p.png.scale,alpha=0.5) # SIZE
		p.png = p.png + geom_text(data=box,aes(color=V4,x=V2,y=max(dm\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.png = p.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.pdf = p.pdf + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,size=0.1*p.pdf.scale,alpha=0.5) # SIZE
		p.pdf = p.pdf + geom_text(data=box,aes(color=V4,x=V2,y=max(dm\$y),label=label),angle=90,vjust=0,hjust=1,size=1) # SIZE
		p.pdf = p.pdf + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		# p2
		p2.png.mod = p2.png.mod + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=1),fill=NA,size=0.5*p.png.scale,alpha=0.5) # SIZE
		p2.png.mod = p2.png.mod + geom_text(data=box,aes(color=V4,x=V2,y=1,label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p2.png.mod = p2.png.mod + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p2.pdf.mod = p2.pdf.mod + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=1),fill=NA,size=0.1*p.pdf.scale,alpha=0.5) # SIZE
		p2.pdf.mod = p2.pdf.mod + geom_text(data=box,aes(color=V4,x=V2,y=1,label=label),angle=90,vjust=0,hjust=1,size=1) # SIZE
		p2.pdf.mod = p2.pdf.mod + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		# p3
		p3.CGprof.col = c(p3.CGprof.col,p.col2)
		p3.CGprof.mod.png = p3.CGprof.mod.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=1),fill=NA,size=0.5*p.png.scale,alpha=0.5) # SIZE
		p3.CGprof.mod.png = p3.CGprof.mod.png + geom_text(data=box,aes(color=V4,x=V2,y=1,label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p3.CGprof.mod.png = p3.CGprof.mod.png + scale_fill_manual(values=c(p3.CGprof.col)) + scale_color_manual(values=c(p3.CGprof.col))

		p3.CGprof.col = c(p3.CGprof.col,p.col2)
		p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=1),fill=NA,size=0.1*p.pdf.scale,alpha=0.5) # SIZE
		p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + geom_text(data=box,aes(color=V4,x=V2,y=1,label=label),angle=90,vjust=0,hjust=1,size=1) # SIZE
		p3.CGprof.mod.pdf = p3.CGprof.mod.pdf + scale_fill_manual(values=c(p3.CGprof.col)) + scale_color_manual(values=c(p3.CGprof.col))

		# p.heatmaponly
		p.heatmaponly2 = p.heatmaponly
		p.heatmaponly2 = p.heatmaponly2 + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm\$y)),fill=NA,size=0.5*p.png.scale,alpha=0.5) # SIZE
		p.heatmaponly2 = p.heatmaponly2 + geom_text(data=box,aes(color=V4,x=V2,y=max(dm\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.heatmaponly2 = p.heatmaponly2 + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.conly.png = p.conly.png + geom_rect(data=box2,aes(color=V4,xmin=beg,xmax=end,ymin=0,ymax=max(dm\$y)),fill=NA,size=0.5*p.png.scale,alpha=0.5)
		p.conly.png = p.conly.png + geom_text(data=box2,aes(color=V4,x=beg,y=max(dm\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.conly.png = p.conly.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		# p.rand.300.png 240717
		p.rand.300.png = p.rand.300.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.300\$y)),fill=NA,size=0.5*p.png.scale) # SIZE
		p.rand.300.png = p.rand.300.png + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.300\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.rand.300.png = p.rand.300.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.heatmaponly2.rand.300.png = p.heatmaponly.rand.300
		p.heatmaponly2.rand.300.png = p.heatmaponly2.rand.300.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.300\$y)),fill=NA,size=0.5*p.png.scale) # SIZE
		p.heatmaponly2.rand.300.png = p.heatmaponly2.rand.300.png + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.300\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.heatmaponly2.rand.300.png = p.heatmaponly2.rand.300.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.rand.300.pdf = p.rand.300.pdf + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.300\$y)),fill=NA,size=0.1*p.pdf.scale) # SIZE
		p.rand.300.pdf = p.rand.300.pdf + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.300\$y),label=label),angle=90,vjust=0,hjust=1,size=1) # SIZE
		p.rand.300.pdf = p.rand.300.pdf + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.heatmaponly2.rand.300.pdf = p.heatmaponly.rand.300
		p.heatmaponly2.rand.300.pdf = p.heatmaponly2.rand.300.pdf + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.300\$y)),fill=NA,size=0.1*p.pdf.scale) # SIZE
		p.heatmaponly2.rand.300.pdf = p.heatmaponly2.rand.300.pdf + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.300\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.heatmaponly2.rand.300.pdf = p.heatmaponly2.rand.300.pdf + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))
	}
	";

	$Rbox_nopk .= "


		print(\"!!! NOPK\")
		print(names(p.col2))
		print(p.col2)
		box2 = box
		box2 = merge(box2,subset(dmc2.rand.1000,select=c(V2,beg)),by=c(\"V2\"))
		box2 = merge(box2,subset(dmc2.rand.1000,select=c(V3,end)),by=c(\"V3\"))
      box2\$x = 0
      box2\$y = 0
      box2\$variable = 0
      box\$x = 0
      box\$y = 0
      box\$variable = 0
		p.rand.1000.png = p.rand.1000.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,size=0.5*p.png.scale) # SIZE
		p.rand.1000.png = p.rand.1000.png + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.1000\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.rand.1000.png = p.rand.1000.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))


		p.heatmaponly2.rand.1000.png = p.heatmaponly.rand.1000
		p.heatmaponly2.rand.1000.png = p.heatmaponly2.rand.1000.png + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,size=0.5*p.png.scale) # SIZE
		p.heatmaponly2.rand.1000.png = p.heatmaponly2.rand.1000.png + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.1000\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.heatmaponly2.rand.1000.png = p.heatmaponly2.rand.1000.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.rand.1000.pdf = p.rand.1000.pdf + geom_rect(data=box,aes(color=V4,xmin=V2,xmax=V3,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,size=0.1*p.pdf.scale) # SIZE
		p.rand.1000.pdf = p.rand.1000.pdf + geom_text(data=box,aes(color=V4,x=V2,y=max(dm.rand.1000\$y),label=label),angle=90,vjust=0,hjust=1,size=1) # SIZE
		p.rand.1000.pdf = p.rand.1000.pdf + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

		p.rand.1000.conly.png = p.rand.1000.conly.png + geom_rect(data=box2,aes(color=V4,xmin=beg,xmax=end,ymin=0,ymax=max(dm.rand.1000\$y)),fill=NA,size=0.5*p.png.scale) # SIZE
		p.rand.1000.conly.png = p.rand.1000.conly.png + geom_text(data=box2,aes(color=V4,x=beg,y=max(dm.rand.1000\$y),label=label),angle=90,vjust=0,hjust=1,size=3) # SIZE
		p.rand.1000.conly.png = p.rand.1000.conly.png + scale_fill_manual(values=c(p.col2)) + scale_color_manual(values=c(p.col2))

	}
	";
	return($Rbox_peak, $Rbox_nopk);
}

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

totalpeak = $totpeak + $totnopk
if (totalpeak == 0) {
	peakperc = 0
} else {
	peakperc = as.integer(1000*$totpeak / ($totpeak + $totnopk)+0.5)/10
}
p.rand.100 = ggplot(dm.rand.100,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_void() +
	 theme(legend.position=\"none\") +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 ggtitle(paste(\"$mytitle\",\"\\n\",peakperc,\"\% $totpeak/\",totalpeak,sep=\"\"))

p.heatmaponly.rand.100 = ggplot(dm.rand.100,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_void() +
	 theme(legend.position=\"none\") +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0))

	p.rand.100.pdf = p.rand.100 + theme(plot.title = element_text(size = 20*p.pdf.scale))
	p.rand.100.png = p.rand.100 + theme(plot.title = element_text(size = 20*p.png.scale))
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
dm.rand.1000\$variable = as.numeric(as.character(dm.rand.1000\$variable))

if (convwant == \"CH\") {
	dmc.rand.1000 = myseq[myseq\$value == \"C\" & myseq\$value2 != \"G\",]
} else if (convwant == \"CG\") {
	dmc.rand.1000 = myseq[myseq\$value == \"C\",]
} else if (convwant == \"GH\") {
	dmc.rand.1000 = myseq[myseq\$value == \"G\" & myseq\$value0 != \"C\",]
} else if (convwant == \"GC\") {
	dmc.rand.1000 = myseq[myseq\$value == \"G\",]
}
dmc.rand.1000 = dmc.rand.1000[order(dmc.rand.1000\$variable),]
dmc.rand.1000\$x= seq(1,dim(dmc.rand.1000)[1])
dmc.rand.1000 = subset(dmc.rand.1000,select=c(variable,x))
dmc.rand.1000 = dmc.rand.1000[order(dmc.rand.1000\$variable),]
i0 = dmc.rand.1000[dim(dmc.rand.1000)[1],,drop=F]
dmc2.rand.1000 = data.frame()
for (i in seq(max(dmc.rand.1000\$variable),1,-1)) {
   if (dim(dmc.rand.1000[dmc.rand.1000\$variable == i,])[1] == 0) {
      icurr = i0
      icurr\$variable = i
      dmc2.rand.1000 = rbind(dmc2.rand.1000,icurr)
   } else {
      i0 = dmc.rand.1000[dmc.rand.1000\$variable == i,,drop=F]
      dmc2.rand.1000 = rbind(dmc2.rand.1000,i0)
   }
}
dmc2.rand.1000 = dmc2.rand.1000[order(dmc2.rand.1000\$variable),]
rownames(dmc2.rand.1000) = seq(1,dim(dmc2.rand.1000)[1])
dmc2.rand.1000\$beg = dmc2.rand.1000\$x
dmc2.rand.1000\$end = dmc2.rand.1000\$x
dmc2.rand.1000\$V2 = dmc2.rand.1000\$variable
dmc2.rand.1000\$V3 = dmc2.rand.1000\$variable

totalpeak = $totpeak + $totnopk
if (totalpeak == 0) {
	peakperc = 0
} else {
	peakperc = as.integer(1000*$totpeak / ($totpeak + $totnopk)+0.5)/10
}
p.rand.1000 = ggplot(dm.rand.1000,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_void() +
	 theme(legend.position=\"none\") +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 ggtitle(paste(\"$mytitle\",\"\\n\",peakperc,\"\% $totpeak/\",totalpeak,sep=\"\"))

p.heatmaponly.rand.1000 = ggplot(dm.rand.1000,aes(variable,y)) +  
	 geom_tile(aes(fill=as.factor(value))) + 
	 theme_void() +
	 theme(legend.position=\"none\") +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0))

dm.rand.1000.conly = dm.rand.1000[dm.rand.1000\$variable \%in\% dmc.rand.1000\$variable,]
dm.rand.1000.conly = merge(dm.rand.1000.conly,dmc.rand.1000,by=\"variable\")

p.rand.1000.conly = ggplot(dm.rand.1000.conly,aes(x,y)) +  
	geom_tile(aes(fill=as.factor(value))) +
	 theme_void() +
	 theme(legend.position=\"none\") +
	 scale_fill_manual(values=c(p.col)) +
	 scale_color_manual(values=c(p.col)) +
	 scale_x_continuous(expand = c(0,0)) + 
	 scale_y_continuous(expand = c(0,0)) +
	 labs(x=NULL,y=NULL)

	p.rand.1000.pdf = p.rand.1000 + theme(plot.title = element_text(size = 20*p.pdf.scale))
	p.rand.1000.png = p.rand.1000 + theme(plot.title = element_text(size = 20*p.png.scale))
	p.rand.1000.conly.png = p.rand.1000.conly

";

# -------------------- $R->{mainplot_nopk}
	$R->{mainplot_nopk} .= "

##########################
# Main Plot Part 1 nopk

dm = melt(df,id.vars=c(\"V1\",\"y\"))
dm\$variable = as.numeric(as.character(dm\$variable))

if (convwant == \"CH\") {
	dmc = myseq[myseq\$value == \"C\" & myseq\$value2 != \"G\",]
} else if (convwant == \"CG\") {
	dmc = myseq[myseq\$value == \"C\",]
} else if (convwant == \"GH\") {
	dmc = myseq[myseq\$value == \"G\" & myseq\$value0 != \"C\",]
} else if (convwant == \"GC\") {
	dmc = myseq[myseq\$value == \"G\",]
}
dmc = dmc[order(dmc\$variable),]
dmc\$x= seq(1,dim(dmc)[1])
dmc = dmc[order(dmc\$variable),]
dmc = subset(dmc,select=c(variable,x))

i0 = dmc[dim(dmc)[1],,drop=F]
dmc2 = data.frame()
for (i in seq(max(dmc\$variable),1,-1)) {
   if (dim(dmc[dmc\$variable == i,])[1] == 0) {
      icurr = i0
      icurr\$variable = i
      dmc2 = rbind(dmc2,icurr)
   } else {
      i0 = dmc[dmc\$variable == i,,drop=F]
      dmc2 = rbind(dmc2,i0)
   }
}
dmc2 = dmc2[order(dmc2\$variable),]
rownames(dmc2) = seq(1,dim(dmc2)[1])
dmc2\$beg = dmc2\$x
dmc2\$end = dmc2\$x
dmc2\$V2 = dmc2\$variable
dmc2\$V3 = dmc2\$variable

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
	p = ggplot(dm.temp,aes(variable,y)) +  
		 geom_tile(aes(fill=as.factor(value))) + 
		 theme_void() +
		 theme(legend.position=\"none\") +
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0))
		 ggtitle(paste(\"PLOT #\",i,\" (total peak=\",$totpeak,\"; total nopk=\",$totnopk,\")\",sep=\"\"))
	plot_list[[i]] = p


	dm.conly = dm.temp[dm.temp\$variable \%in\% dmc\$variable,]
	dm.conly = merge(dm.conly,dmc,by=\"variable\")
	p.conly = ggplot(dm.conly,aes(x,y)) +  
		geom_tile(aes(fill=as.factor(value))) +
		 theme_void() +
		 theme(legend.position=\"none\") +
		 scale_fill_manual(values=c(p.col)) +
		 scale_color_manual(values=c(p.col)) +
		 scale_x_continuous(expand = c(0,0)) + 
		 scale_y_continuous(expand = c(0,0)) +
		 labs(x=NULL,y=NULL)
	p.png = p
	p.pdf = p
	p.heatmaponly = p
	p.conly.png = p.conly
}
";


# -------------------- $R->{mainplotClusterAddition}
	$R->{mainplotClusterAddition} = "

# Main Plot Cluster Addition
# geom rect default size = 0.5
# geom text default size = 10
p.png = p.png + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.png.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=paste(clust-10,\"(\",id2,\")\",sep=\"\")),hjust=0,size=5*p.png.scale) +
	 theme(plot.title = element_text(size = 20*p.png.scale))

p.pdf = p.pdf + 
	 geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,size=0.5*p.pdf.scale) + # SIZE
	 geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=paste(clust-10,\"(\",id2,\")\",sep=\"\")),hjust=0,size=5*p.pdf.scale) +
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
if (length(df2[df2 < 4]) > 0) {
   df2[df2 < 4] = 0;
}
if (length(df2[df2 >= 4 & df2 < 6]) > 0) {
   df2[df2 >= 4 & df2 < 6] = 1
}
if (length(df2[df2 >= 6 & df2 < 8]) > 0) {
   df2[df2 >= 6 & df2 < 8] = 2
}
if (length(df2[df2 >= 8 & df2 <= 9]) > 0) {
   df2[df2 >= 8 & df2 <= 9] = 3
}
df2temp = data.frame(x=seq(1,dim(df2)[2]),C=apply(df2,2,max))
df2temp2 = df2
if (length(df2temp2[df2temp2 == 1]) > 0) {
   df2temp2[df2temp2 == 1] = 0;
}
if (length(df2temp2[df2temp2 >= 2]) > 0) {
   df2temp2[df2temp2 >= 2] = 1;
}
df2temp\$meanall = apply(df2temp2,2,mean)

df2temp2 = df2
if (length(df2temp2[df2temp2 == 1]) > 0) {
   df2temp2[df2temp2 == 1] = 0;
}
if (length(df2temp2[df2temp2 >= 3]) > 0) {
   df2temp2[df2temp2 >= 3] = 1;
}
df2temp\$meanpeak = apply(df2temp2,2,mean)

if (convwant == \"CG\") {
	df2temp = df2temp[df2temp\$x \%in\% myseq[myseq\$type == \"CG\" | myseq\$type == \"CH\",]\$variable,]
} else if (convwant == \"GC\") {
	df2temp = df2temp[df2temp\$x \%in\% myseq[myseq\$type == \"GC\" | myseq\$type == \"GH\",]\$variable,]
} else if (convwant == \"CH\") {
	df2temp = df2temp[df2temp\$x \%in\% myseq[myseq\$type == \"CH\",]\$variable,]
} else if (convwant == \"GH\") {
	df2temp = df2temp[df2temp\$x \%in\% myseq[myseq\$type == \"GH\",]\$variable,]
}

# using meanall
if (dim(df2temp[df2temp\$meanall > 0,])[1] > 0) {
   df2 = data.frame(x=df2temp\$x,y=df2temp\$meanall,x2=df2temp\$x,y2=df2temp\$meanall)
} else {
   df2 = data.frame(x=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y=0, x2=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y2=0);
}
#old
write.table(df2,file=\"$currFile.fixedwig\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
write.table(df2temp,file=\"$currFile.fixedwig2\",quote=F,row.names=F,col.names=T,sep=\"\\t\")

# P2 % Conversion XY Plot
	p2.png.scale = p.png.scale
	p2.pdf.scale = p.pdf.scale

p2.png = 
	ggplot(df2,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=0.5*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),lwd=0.5) + # SIZE
	theme_void() +
	theme(legend.position=\"none\") +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	coord_cartesian(ylim=c(-0.05,1.05)) +
	annotate(geom='segment',x=0.01*max(df2\$x2),xend=0.01*max(df2\$x2),y=0,yend=1,color=\"black\")
	for (i in seq(0,1,0.25)) {
		p2.png = p2.png + 
			annotate(geom='segment',x=0.01*max(df2\$x2),xend=0.02*max(df2\$x2),y=i,yend=i,color=\"black\") +
			annotate(geom='text',x=0.02*max(df2\$x2),y=i,label=paste(i*100,\" \%\",sep=\"\"),size=5*p2.png.scale,hjust=0)
	}

p2.png.mod = p2.png

p2.pdf = 
	ggplot(df2,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),pch=\".\",size=0.1) +  # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),lwd=0.1) + # SIZE 
	theme_void() +
	theme(legend.position=\"none\") +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	coord_cartesian(ylim=c(-0.05,1.05)) +
	annotate(geom='segment',x=0.01*max(df2\$x2),xend=0.01*max(df2\$x2),y=0,yend=1,color=\"black\")
	for (i in seq(0,1,0.25)) {
		p2.pdf = p2.pdf + 
			annotate(geom='segment',x=0.01*max(df2\$x2),xend=0.02*max(df2\$x2),y=i,yend=i,color=\"black\") +
			annotate(geom='text',x=0.02*max(df2\$x2),y=i,label=paste(i*100,\" \%\",sep=\"\"),hjust=0,size=1)
	}

p2.pdf.mod = p2.pdf
";

	$R->{secondplotConversionGraph_rand_100} = "

# Calculate % Conversion
df2.rand.100 = subset(df.rand.100,select=c(-V1,-y));
if (length(df2.rand.100[df2.rand.100 < 4]) > 0) {
   df2.rand.100[df2.rand.100 < 4] = 0;
}
if (length(df2.rand.100[df2.rand.100 >= 4 & df2.rand.100 < 6]) > 0) {
   df2.rand.100[df2.rand.100 >= 4 & df2.rand.100 < 6] = 1
}
if (length(df2.rand.100[df2.rand.100 >= 6 & df2.rand.100 < 8]) > 0) {
   df2.rand.100[df2.rand.100 >= 6 & df2.rand.100 < 8] = 2
}
if (length(df2.rand.100[df2.rand.100 >= 8 & df2.rand.100 <= 9]) > 0) {
   df2.rand.100[df2.rand.100 >= 8 & df2.rand.100 <= 9] = 3
}
df2.rand.100temp = data.frame(x=seq(1,dim(df2.rand.100)[2]),C=apply(df2.rand.100,2,max))
df2.rand.100temp2 = df2.rand.100
if (length(df2.rand.100temp2[df2.rand.100temp2 == 1]) > 0) {
   df2.rand.100temp2[df2.rand.100temp2 == 1] = 0;
}
if (length(df2.rand.100temp2[df2.rand.100temp2 >= 2]) > 0) {
   df2.rand.100temp2[df2.rand.100temp2 >= 2] = 1;
}
df2.rand.100temp\$meanall = apply(df2.rand.100temp2,2,mean)

df2.rand.100temp2 = df2.rand.100
if (length(df2.rand.100temp2[df2.rand.100temp2 == 1]) > 0) {
   df2.rand.100temp2[df2.rand.100temp2 == 1] = 0;
}
if (length(df2.rand.100temp2[df2.rand.100temp2 >= 3]) > 0) {
   df2.rand.100temp2[df2.rand.100temp2 >= 3] = 1;
}
df2.rand.100temp\$meanpeak = apply(df2.rand.100temp2,2,mean)

# using meanall
df2.rand.100temp = df2.rand.100temp[df2.rand.100temp\$C != 0,]
if (dim(df2.rand.100temp[df2.rand.100temp\$meanall > 0,])[1] > 4) {
   df2.rand.100 = data.frame(x=df2.rand.100temp\$x,y=df2.rand.100temp\$meanall,x2=df2.rand.100temp\$x,y2=df2.rand.100temp\$meanall)
} else {
   df2.rand.100 = data.frame(x=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y=0, x2=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y2=0);
}

# P2 % Conversion XY Plot
p2.rand.png.scale = p.png.scale
p2.rand.pdf.scale = p.pdf.scale

p2.rand.100.png = 
	ggplot(df2.rand.100,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_void() +
	theme(legend.position=\"none\") +
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
	theme_void() +
	theme(legend.position=\"none\") +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=1,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{secondplotConversionGraph_rand_1000} = "

# Calculate % Conversion
df2.rand.1000 = subset(df.rand.1000,select=c(-V1,-y));
if (length(df2.rand.1000[df2.rand.1000 < 4]) > 0) {
   df2.rand.1000[df2.rand.1000 < 4] = 0;
}
if (length(df2.rand.1000[df2.rand.1000 >= 4 & df2.rand.1000 < 6]) > 0) {
   df2.rand.1000[df2.rand.1000 >= 4 & df2.rand.1000 < 6] = 1
}
if (length(df2.rand.1000[df2.rand.1000 >= 6 & df2.rand.1000 < 8]) > 0) {
   df2.rand.1000[df2.rand.1000 >= 6 & df2.rand.1000 < 8] = 2
}
if (length(df2.rand.1000[df2.rand.1000 >= 8 & df2.rand.1000 <= 9]) > 0) {
   df2.rand.1000[df2.rand.1000 >= 8 & df2.rand.1000 <= 9] = 3
}
df2.rand.1000temp = data.frame(x=seq(1,dim(df2.rand.1000)[2]),C=apply(df2.rand.1000,2,max))
df2.rand.1000temp2 = df2.rand.1000
if (length(df2.rand.1000temp2[df2.rand.1000temp2 == 1]) > 0) {
   df2.rand.1000temp2[df2.rand.1000temp2 == 1] = 0;
}
if (length(df2.rand.1000temp2[df2.rand.1000temp2 >= 2]) > 0) {
   df2.rand.1000temp2[df2.rand.1000temp2 >= 2] = 1;
}
df2.rand.1000temp\$meanall = apply(df2.rand.1000temp2,2,mean)

df2.rand.1000temp2 = df2.rand.1000
if (length(df2.rand.1000temp2[df2.rand.1000temp2 == 1]) > 0) {
   df2.rand.1000temp2[df2.rand.1000temp2 == 1] = 0;
}
if (length(df2.rand.1000temp2[df2.rand.1000temp2 >= 3]) > 0) {
   df2.rand.1000temp2[df2.rand.1000temp2 >= 3] = 1;
}
df2.rand.1000temp\$meanpeak = apply(df2.rand.1000temp2,2,mean)

# using meanall
df2.rand.1000temp = df2.rand.1000temp[df2.rand.1000temp\$C != 0,]
if (dim(df2.rand.1000temp[df2.rand.1000temp\$meanall > 0,])[1] > 4) {
   df2.rand.1000 = data.frame(x=df2.rand.1000temp\$x,y=df2.rand.1000temp\$meanall,x2=df2.rand.1000temp\$x,y2=df2.rand.1000temp\$meanall)
} else {
   df2.rand.1000 = data.frame(x=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y=0, x2=seq(1,dim(subset(df,select=c(-V1,-y)))[2]), y2=0);
}
 
# P2 % Conversion XY Plot
p2.rand.png.scale = p.png.scale
p2.rand.pdf.scale = p.pdf.scale

p2.rand.1000.png = 
	ggplot(df2.rand.1000,aes(x2,y2)) + 
	geom_point(aes(x=x,y=y),size=1*p2.png.scale) + # SIZE
	geom_line(color=rgb(1,0,0,alpha=1),size=1*p2.png.scale) + # SIZE
	theme_void() +
	theme(legend.position=\"none\") +
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
	theme_void() +
	theme(legend.position=\"none\") +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
	annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0.75,label=\"-  75 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=5,label=\"-  50 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0.25,label=\"-  25 \%\",size=1,hjust=0) +
	annotate(geom='text',x=10,y=0,label=\"-   0\%\",size=1,hjust=0) +
	coord_cartesian(ylim=c(-0.05,1.05))

";

	$R->{Scale} = "

ratio1 = max(30,dim(df)[1]) / (max(30,dim(df)[1]) + 31.25 + 26.5625)
ratio2 = 31.25      / (max(30,dim(df)[1]) + 31.25 + 26.5625)
ratio3 = 26.5625    / (max(30,dim(df)[1]) + 31.25 + 26.5625)
mynrow = 2
totalheight = max(30,dim(df)[1]) + 31.25
totalheight2 = max(30,dim(df)[1]) + 31.25 + 26.5625
totalwidth  = dim(df)[2] * 0.125
totalratio  = c(ratio1, ratio2)
totalheight = totalheight + 26.5625
totalratio = c(totalratio, ratio3)
mynrow = 3

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



############ Declare ###########
$R->{PNG} = "";
$R->{PNG_peak_300} = "";
$R->{PNG_nopk} = "";
$R->{PNG_nopk_rand_100} = "";
$R->{PNG_nopk_ALL} = "";

$R->{PDF} = "";
$R->{PDF_peak_300} = "";
$R->{PDF_nopk} = "";
$R->{PDF_nopk_rand_100} = "";
$R->{PDF_nopk_ALL} = "";
################################

$R->{PNG_peak_300} = "

	# PNG HEATMAP ONLY SUBSAMPLE
	totalheight = dim(df.rand.300)[1] * myscale
	pngout_heatmap_only = \"$pngoutFolder/HEATMAP300/$pngoutFilename.heatmap.300.png\"
	print(paste(\"Making\",pngout_heatmap_only))
	png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
	print(p.heatmaponly.rand.300)
	dev.off()
	
	pngout_heatmaponly2 = \"$pngoutFolder/HEATMAPBOX300/$pngoutFilename.heatmapbox.300.png\"
	print(paste(\"Making\",pngout_heatmaponly2))
	png(pngout_heatmaponly2,width=min(30000,totalwidth),height=min(30000,totalheight))
	print(p.heatmaponly2.rand.300.png)
	dev.off()

";

$R->{PNG} = "

	# PNG
	print(\"$pngout\")
	png(\"$pngout\",width=min(30000,totalwidth),height=min(30000,totalheight))
	grid.arrange(p.png,p2.png.mod,p3.CGprof.mod.png,ncol=1,nrow=3,heights=totalratio)
	dev.off()
	
	# PNG HEATMAP ONLY
	totalheight = max(30,dim(df)[1]) * myscale
	pngout_heatmap_only = \"$pngoutFolder/HEATMAP/$pngoutFilename.heatmap.png\"
	print(pngout_heatmap_only)
	png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
	print(p.heatmaponly)
	dev.off()
	
	# PNG HEATMAP2 (+ box)
	pngout_heatmaponly2 = \"$pngoutFolder/HEATMAP2/$pngoutFilename.heatmapbox.png\"
	print(pngout_heatmaponly2)
	png(pngout_heatmaponly2,width=min(30000,totalwidth),height=min(30000,totalheight))
	grid.arrange(p.heatmaponly2)
	dev.off()
	
	# PNG all Conv
	pngout_peak_all_c_conv = \"$pngoutFolder/CONV/$pngoutFilename.c_conv.png\"
	print(pngout_peak_all_c_conv)
	png(pngout_peak_all_c_conv,width=min(30000,totalwidth),height=31.25*myscale)
	grid.arrange(p2.png)
	dev.off()
	
	# PNG all Conv
	pngout_peak_cgprof = \"$pngoutFolder/CGPROF/$pngoutFilename.cgprof.png\"
	print(pngout_peak_cgprof)
	png(pngout_peak_cgprof,width=min(30000,totalwidth),height=31.25*myscale)
	grid.arrange(p3.CGprof)
	dev.off()
	
	totalwidth.conly = max(dm.conly\$x)*4
	totalheight.conly = max(dm.conly\$y)*4
	# PNG Conly
	pngout_peak_all_conly = \"$pngoutFolder/CONLY/$pngoutFilename.conly.png\"
	print(pngout_peak_all_conly)
	png(pngout_peak_all_conly,width=min(30000,totalwidth.conly),height=min(30000,totalheight.conly))
	grid.arrange(p.conly.png)
	dev.off()

";

$R->{PNG_nopk} = "

	# PNG
	totalheight = (dim(df.rand.1000)[1] + 31.25) * myscale
	print(\"$pngout\")
	png(\"$pngout\",width=min(30000,totalwidth),height=min(30000,totalheight))
	grid.arrange(p.rand.1000.png,p2.rand.1000.png,ncol=1,nrow=mynrow,heights=totalratio)
	dev.off()
	
	# PNG HEATMAP ONLY
	totalheight = dim(df.rand.1000)[1] * myscale
	pngout_heatmap_only = \"$pngoutFolder/HEATMAP/$pngoutFilename.heatmap.png\"
	print(pngout_heatmap_only)
	png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
	print(p.heatmaponly.rand.1000)
	dev.off()
	
	# PNG all Conv
	pngout_nopk_all_c_conv = \"$pngoutFolder/CONV/$pngoutFilename.c_conv.png\"
	print(pngout_nopk_all_c_conv)
	png(pngout_nopk_all_c_conv,width=min(30000,totalwidth),height=31.25*myscale)
	grid.arrange(p2.png)
	dev.off()
	
	totalwidth.conly = max(dm.rand.1000.conly\$x)*4
	totalheight.conly = max(dm.rand.1000.conly\$y)*4
	# PNG all Conv
	pngout_nopk_all_conly = \"$pngoutFolder/CONLY/$pngoutFilename.conly.png\"
	print(pngout_nopk_all_conly)
	png(pngout_nopk_all_conly,width=min(30000,totalwidth.conly),height=min(30000,totalheight.conly))
	grid.arrange(p.rand.1000.conly.png)
	dev.off()

";

$R->{PNG_nopk_rand_100} = "

	pngout_nopk_rand_100 = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.png\"
	# PNG
	totalheight = (dim(df.rand.100)[1] + 31.25) * myscale
	png(pngout_nopk_rand_100,width=min(30000,totalwidth),height=min(30000,totalheight))
	grid.arrange(p.rand.100.png,p2.rand.100.png,ncol=1,nrow=mynrow,heights=totalratio)
	dev.off()
	
	# PNG HEATMAP ONLY
	totalheight = dim(df.rand.100)[1] * myscale
	pngout_heatmap_only = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.heatmap.png\"
	png(pngout_heatmap_only,width=min(30000,totalwidth),height=min(30000,totalheight))
	print(p.heatmaponly.rand.100)
	dev.off()
	
	# PNG all Conv
	pngout_nopk_all_c_conv = \"$pngoutFolder/ALL/$pngoutFilename.RAND.100.c_conv.png\"
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
			png(pngout_nopk,width=min(30000,totalwidth),height=min(30000,currtotalheight))
			grid.arrange(plot_list[[i]],p2,ncol=1,nrow=mynrow,heights=currtotalratio)
			dev.off()
		} else {
			currtotalheight = totalheight_nopk
			png(pngout_nopk,width=min(30000,totalwidth),height=min(30000,currtotalheight))
			grid.arrange(plot_list[[i]],ncol=1)
			dev.off()
		}
	}
";


####### PDF ########

my $PDFSCALE = 200;

$R->{PDF_peak_300} = "

	# PDF HEATMAP ONLY SUBSAMPLE 240717
	currheight = dim(df.rand.300)[1] * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_heatmap_only = \"$pdfoutFolder/HEATMAP300/$pdfoutFilename.heatmap.300.pdf\"
	print(paste(\"Making\",pdfout_heatmap_only))
	print(paste(\"currheight\",currheight,currwidth))
	pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
	print(p.heatmaponly.rand.300)
	dev.off()
	
	# PDF HEATMAP2 ONLY SUBSAMPLE 240717
	currheight = dim(df.rand.300)[1] * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_heatmaponly2 = \"$pdfoutFolder/HEATMAPBOX300/$pdfoutFilename.heatmapbox.300.pdf\"
	print(paste(\"Making\",pdfout_heatmaponly2))
	print(paste(\"currheight\",currheight,currwidth))
	pdf(pdfout_heatmaponly2,width=min(30000,currwidth),height=min(30000,currheight))
	print(p.heatmaponly2.rand.300.pdf)
	dev.off()
";

$R->{PDF} = "

	# PDF
	currheight = totalheight / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	print(\"$pdfout\")
	pdf(\"$pdfout\",width=min(30000,currwidth),height=min(30000,currheight))
	grid.arrange(p.pdf,p2.pdf.mod,p3.CGprof.mod.pdf,ncol=1,nrow=3,heights=totalratio)
	dev.off()
	
	# PDF all Conv
	pdfout_peak_cgprof = \"$pdfoutFolder/CGPROF/$pdfoutFilename.cgprof.pdf\"
	print(pdfout_peak_cgprof)
	pdf(pdfout_peak_cgprof,width=min(30000,currwidth),height=min(30000,currheight))
	grid.arrange(p3.CGprof.mod.pdf)
	dev.off()
	
	# PDF HEATMAP ONLY
	currheight = dim(df)[1] * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_heatmap_only = \"$pdfoutFolder/HEATMAP/$pdfoutFilename.heatmap.pdf\"
	pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
	print(p.heatmaponly)
	dev.off()
	
	# PDF all Conv
	currheight = 31.25 * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_peak_all_c_conv = \"$pdfoutFolder/CONV/$pdfoutFilename.c_conv.pdf\"
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
	pdfout_heatmap_only = \"$pdfoutFolder/HEATMAP/$pdfoutFilename.heatmap.pdf\"
	pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
	print(p.heatmaponly.rand.1000)
	dev.off()
	
	# PDF all Conv
	currheight = 31.25 * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_nopk_all_c_conv = \"$pdfoutFolder/CONV/$pdfoutFilename.c_conv.pdf\"
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
	pdfout_heatmap_only = \"$pdfoutFolder/HEATMAP/$pdfoutFilename.RAND.100.heatmap.pdf\"
	pdf(pdfout_heatmap_only,width=min(30000,currwidth),height=min(30000,currheight))
	print(p.heatmaponly.rand.100)
	dev.off()
	
	
	# PDF all Conv
	currheight = 31.25 * myscale / $PDFSCALE
	currwidth = totalwidth / $PDFSCALE
	pdfout_nopk_all_c_conv = \"$pdfoutFolder/CONV/$pdfoutFilename.RAND.100.c_conv.pdf\"
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
