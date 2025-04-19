#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_w $opt_g $opt_G $opt_v $opt_n $opt_r $opt_R $opt_B $opt_c $opt_F $opt_0 $opt_J $opt_C $opt_f $opt_p);
getopts("n:vg:w:G:r:R:B:cfF0J:Cp");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
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

my $footPeakGraphMainScript = "$footLoop_script_folder/bin/footPeak_graph_main.pl";
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
	my %newfiles;
	my ($footPeak_logFile) = "$resDir/footPeak_logFile.txt";
	open (my $outLog, ">", "$resDir/footPeak_graph_logFile.txt") or die "\n\nFailed to write to $resDir/footPeak_graph_logFile.txt: $!\n\n";

	LOG($outLog, ">footPeak_graph.pl version $version\n");
	LOG($outLog, ">UUID: $uuid\n", "NA");
	LOG($outLog, ">Date: $date\n", "NA");
	LOG($outLog, ">Run script: $0 -n $opt_n\n", "NA");
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
	my @newFiles;
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
				next if (defined $opt_G and $peakFile !~ /$opt_G/);
				my $peakflag = getFlag($peakFile, "Pos", $strand, $type);
				my $nopkflag = getFlag($nopkFile, "Pos", $strand, $type);
				if ($peakflag !~ /RCONV/) {
					LOG($outLog, "${LPR}PNGOUT $peakflag:$N\n");
					LOG($outLog, "scp mitochi\@franklin.hpc.ucdavis.edu:/$resDir2/PNG/$peakflag/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.PEAK.out.$peakflag.png.ALL.conly.png ./\n");
					LOG($outLog, "scp mitochi\@franklin.hpc.ucdavis.edu:/$resDir2/PNG/$nopkflag/CONLY/$label\_gene$mygene\_$strand\_$window\_$thres\_$type.NOPK.out.$nopkflag.png.ALL.conly.png ./\n");
				}
				$peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$strand\_$window\_$thres2\_$type.PEAK" if not -e $peakFile;
				
				print "Can't find peakFile $LCY$peakFile$N\n" if not -e $peakFile;
				LOG($outLog, "label=$label, gene$mygene, strand=$strand, peak=$peakFile\n","NA");
				$newfiles{$peakFile} = $mygene;

			   my $RDMADEPNG = $opt_r;
			   $RDMADEPNG = 0 if $peakflag =~ /(RCONV)/ and $opt_r < 2;
			   $RDMADEPNG = 0 if $peakflag =~ /(ALL)/ and $opt_r < 3;
			   $RDMADEPNG = 0 if $peakflag =~ /_C$/ and not defined $opt_c and not defined $opt_C;
			   $RDMADEPNG = 0 if $peakflag !~ /_C$/ and defined $opt_c and not defined $opt_C;

			   my $RDMADEPDF = $opt_R;
			   $RDMADEPDF = 0 if $peakflag =~ /(RCONV)/ and $opt_R < 2;
			   $RDMADEPDF = 0 if $peakflag =~ /(ALL)/ and $opt_R < 3;
			   $RDMADEPDF = 0 if $peakflag =~ /_C$/ and not defined $opt_c and not defined $opt_C;
			   $RDMADEPDF = 0 if $peakflag !~ /_C$/ and defined $opt_c and not defined $opt_C;
				next if ($RDMADEPNG == 0 and $RDMADEPDF == 0);
				push(@newFiles, $peakFile);
			}
		}
	}
	my %Rscripts; 
	my $lastfile = -1; 
	DIELOG($outLog, "\n\nERROR: There is no files in $LCY\%newfiles$N defined!\n") if (keys %newfiles) == 0;
	my $fileCount = 0;
	my $totalFile = (keys %newfiles);
	my $lastGENE = -1;

	my $sbatch_cmd = "$footPeakGraphMainScript -n $resDir -i FILENAME -r $opt_r -R $opt_R";
	$sbatch_cmd .= " -c" if defined $opt_c;
	$sbatch_cmd .= " -C" if defined $opt_C;
	$sbatch_cmd .= " -0" if defined $opt_0;
	$sbatch_cmd .= " -f" if defined $opt_f;
	$sbatch_cmd .= " -F" if defined $opt_F;
	$sbatch_cmd .= " -p" if defined $opt_p;
	$sbatch_cmd .= " -B $boxFile" if defined $opt_B;
	my $force_sbatch = 1 if defined $opt_F or defined $opt_f;
	my $outsbatchDir = "$resDir/.footPeak_graph_sbatch/";
	system("mkdir -p $outsbatchDir") if not -d $outsbatchDir;
	my $debug; my $mem = 4000;
	$debug = 0 if defined $opt_0;

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
		footLoop_sbatch_main($sbatch_cmd, "footPeak_graph", \@newFiles, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir, $mem, $debug);
	}
	else {
	   LOG($outLog, "#Single run\n\n");
	   for (my $i = 0; $i < @newFiles; $i++) {
	      my $newFile = $newFiles[$i];
	      my ($newFilename) = getFilename($newFile, "full");
	      my $newFileDone = "$outsbatchDir/$newFilename.done";
	      my $currcmd = $sbatch_cmd;
	      my $indice = $i + 1;
	      $currcmd =~ s/FILENAME/$newFile/g;
	      $currcmd =~ s/FNINDICE/$indice/g;
	      if (defined $opt_0) {
	         LOG($outLog, date() . "$YW$indice/$totalFile$N $LCY$newFilename$N\n");
	         LOG($outLog, "\t$LGN$currcmd$N\n");
	      }
	      elsif (not -e $newFileDone or defined $opt_F or defined $opt_f) {
	         system($currcmd) == 0 or DIELOG($outLog, "Failed to run cmd: $!\n\n$LCY$currcmd$N\n\n");
	         system("touch $newFileDone") if not -e $newFileDone;
	      }
	      else {
	         LOG($outLog, date() . "${LPR}footPeak.pl$N: $YW$indice/$totalFile$N $LCY$newFilename$N: using previously made peaks\n");
	         LOG($outLog, "Done=$LCY$newFileDone$N)\n","NA");
	      }
	   }
	}
}
