#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_t $opt_R $opt_D $opt_0 $opt_J $opt_F $opt_f $opt_p);
getopts("vd:n:G:t:RD:0J:Ffp");

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

my $date = getDate();
my $uuid = getuuid();

my ($dist, $footPeakFolder) = ($opt_d, $opt_n);
my $toggleRstrand = defined $opt_R ? "Yes" : "No";

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N $LGN-g gene$N $CY-n <footPeak's output folder (footPeak's -o)>$N

-f/F : force rerun
-p   : Parallel run if slurm system exist
-D   : tightness of read placement in each cluster (default: 200)
-R   : toggle to order reads reversely (good foor in vitro REVERSE_ genes)

";

die $usage  unless defined $opt_n and -d $opt_n;

my $outDir = "$footPeakFolder/FOOTCLUST/";
$outDir= getFullpath($outDir);
my $footClustMainScript = "$footLoop_script_folder/bin/footClust_main.pl";
my $windowDiv = defined $opt_D ? $opt_D : 200;
die "\nWindowDiv must be integer!\n" if $windowDiv !~ /^\d+$/ or $windowDiv <= 0;

my $clustThreshold = defined $opt_t ? $opt_t : 3;
die "\n\nERROR! -t must be positive number (pls don't make it scientific format, use a decimal)!\n" if ($clustThreshold !~ /^\d+\.?\d*$/);
makedir($outDir) if not -d $outDir;
die "Failed to create output directory $LCY$outDir$N!\n" unless -d $outDir;
makedir("$outDir/.TEMP");
die "Failed to create output directory $LCY$outDir/.TEMP/$N!\n" unless -d "$outDir/.TEMP";
makedir("$outDir/PNG/");
die "Failed to create output directory $LCY$outDir/PNG$N!\n" unless -d "$outDir/PNG/";

# establish log file
open (my $outLog, ">", "$outDir/logFile_footClust.txt") or die "Failed to create outLog file $outDir/logFile_footClust.txt: $!\n";

my ($max_parallel_run) = defined $opt_J ? $opt_J : 1;

my $footClust_sbatch_logFileDir = "$footPeakFolder/.footClust_sbatch/";
my $footClust_sbatch_2_logFileDir = "$footPeakFolder/.footClust_sbatch_2/";
makedir("$footClust_sbatch_logFileDir") if not -d "$footClust_sbatch_logFileDir";
makedir("$footClust_sbatch_2_logFileDir") if not -d "$footClust_sbatch_2_logFileDir";
LOG($outLog, date() . "Created $LCY$footClust_sbatch_logFileDir$N\n");
LOG($outLog, date() . "Created $LCY$footClust_sbatch_2_logFileDir$N\n");
makedir("$outDir/CLUST_LOCAL") if not -d "$outDir/CLUST_LOCAL";
makedir("$outDir/CLUST_GENOME") if not -d "$outDir/CLUST_GENOME";


# parse footLoop log file
($footPeakFolder) = getFullpath($footPeakFolder);
my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
my ($footLoopFolder);
my ($geneIndexFile);
open (my $infootPeak_logFile, "<", $footPeak_logFile) or DIELOG($outLog, "Failed to read from $footPeak_logFile: $!\n");
while (my $line = <$infootPeak_logFile>) {
	chomp($line);
		#>Options from footLoop.pl logfile=PCB9/0_Fastq/170804_pcb09_BSCset2_ccs_3minFP.fastq.gz.rmdup.fq.gz_MAP85pBUF100/logFile.txt	
	if ($line =~ /^Run script full/) {
		my @runscriptfull = split(" ", $line);
		for (my $i = 0; $i < @runscriptfull; $i++) {
			if ($runscriptfull[$i] eq "-n") {
				($footLoopFolder) = $runscriptfull[$i+1];
				last;

			}
		}
		last;
	}
}
my $footLoopFolderPrint = defined $footLoopFolder ? $footLoopFolder : "UNDEFINED";
die "Cannot find footLoop Folder $LCY$footLoopFolderPrint$N from logfile $footPeak_logFile\n" unless defined $footLoopFolder and -d $footLoopFolder;
close $infootPeak_logFile;
my $footLoop_logFile = $footLoopFolder . "/logFile.txt";
($geneIndexFile) = parse_footLoop_logFile($footLoop_logFile, $date, $uuid, $footLoopFolder);
my %coor;
open (my $inGeneIndexFile, "<", $geneIndexFile) or die;
while (my $line = <$inGeneIndexFile>) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	$coor{uc($gene)}{chr} = $chr;
	$coor{uc($gene)}{beg} = $beg;
	$coor{uc($gene)}{end} = $end;
	$coor{uc($gene)}{strand} = $strand;
}
close $inGeneIndexFile;
LOG($outLog, "cluster threshold = $clustThreshold\ngene Index file = $geneIndexFile\n");

# sanity check -d distance
$dist = 250 if not defined $opt_d;
LOG($outLog, date() . "\n$LRD Error$N: Distance must be digits!\n") and die unless $dist =~ /^\d+$/;


# get .fa file from footPeakFolder and copy
my ($faFile) = <$footPeakFolder/*.fa>;
LOG($outLog, date() . "\n$LRD Error$N: Cannot find any fasta file in $LCY$footPeakFolder$N\n\n") if not defined $faFile or not -e $faFile;
$faFile =~ s/\/+/\//g;
system("/bin/cp $faFile $outDir") == 0 or LOG($outLog, date() . "Failed to copy$LGN $faFile$N to $LCY$outDir$N: $!\n") and die;


# get .local.bed peak files used for clustering
my @local_peak_files_temp = <$footPeakFolder/PEAKS_LOCAL/*PEAK.local.bed>;
makedir("$footPeakFolder/PEAKS_LOCAL/.delete/") if not -d "$footPeakFolder/PEAKS_LOCAL/.delete/";
my @local_peak_files;
LOG($outLog, date() . "Checking PEAKS_LOCAL directory for empty peaks!\n");
my $empty = 0;
my $notempty = 0;
for (my $i = 0; $i < @local_peak_files_temp; $i++) {
	print "$i\t$local_peak_files_temp[$i]\n" if $i % 1e3 == 0;
	my $size = -s $local_peak_files_temp[$i];
	push(@local_peak_files, $local_peak_files_temp[$i]) if $size > 0;
	if ($size == 0) {
		$empty ++;
		system("/bin/mv $local_peak_files_temp[$i] $footPeakFolder/PEAKS_LOCAL/.delete/");
	}
	else {
		$notempty ++;
	}
}
my ($deleted) = `ls $footPeakFolder/PEAKS_LOCAL/.delete/|wc -l ` =~ /^(\d+)/;
LOG($outLog, date() . "Found $LCY$notempty$N good peaks, $LGN$empty$N empty peaks and $LGN$deleted$N deleted peaks\n");
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


# Parse fasta file
my %genes;
open (my $faIn2, "<", $faFile) or die;
open (my $faOut, ">", "$faFile.footClustfastafile") or die;
my $fasta2 = new FAlite($faIn2);
while (my $entry = $fasta2->nextEntry()) {
	my $def = uc($entry->def);
	my $seq = $entry->seq;
	print $faOut "$def\n$seq\n";
	my $length = length($seq);
	$def =~ s/>//;
	$genes{uc($def)} = $length;
}
close $faIn2;
close $faOut;
system("/bin/rm $faFile.footClustfastafile.fai") if -e "$faFile.footClustfastafile.fai";
system("samtools faidx $faFile.footClustfastafile") == 0 or LOG($outLog, "Failed to samtools faidx $faFile.footClustfastafile: $!\n") and exit 0;

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
my $curr_gene = -1;
my $type_count = -1;
my @Rscript;
my $totalpeakFiles = @local_peak_files;
my @neworigFiles;
for (my $i = 0; $i < @local_peak_files; $i++) {
	my $local_peak_file = $local_peak_files[$i];
	my $peakFile = $local_peak_file;
	if (defined $opt_G) {
		if ($local_peak_file !~ /$opt_G/i) {
			next;
		}
		else {
			LOG($outLog, date() . "${LPR}footClust.pl$N: $LGN$i/$totalpeakFiles$N ${LGN}Processing$N $LCY$local_peak_file$N as it contain $LGN-G $opt_G$N\n");
		}
	}
	else {
		LOG($outLog, date() . "${LPR}footClust.pl$N: $LGN$i/$totalpeakFiles$N ${LGN}Processing$N $LCY$local_peak_file$N\n");
	}

	push(@neworigFiles, $peakFile);
}
my $totalorigFile = @neworigFiles;

my @cmd = ("FILENAME", $faFile, $label, $footLoopFolder, $geneIndexFile, $totalpeakFiles);
my $cmd = "$footClustMainScript";
$cmd .= " -n $footPeakFolder";
my $cmd2 = "-J $max_parallel_run -t $clustThreshold -D $windowDiv -d $dist";
$cmd2 .= " -R" if defined $opt_R;
$cmd2 .= " -F" if defined $opt_F;
$cmd2 .= " -f" if defined $opt_f;
$cmd2 .= " -G $opt_G" if defined $opt_G;
my $runScript = $cmd2;
$cmd .= " $cmd2 " . join(" ", @cmd);
my $sbatch_cmd = $cmd; 
my $force_sbatch = 1 if defined $opt_F;
my $outsbatchDir = "$footPeakFolder/.footClust_sbatch/";
system("mkdir -p $outsbatchDir") if not -d $outsbatchDir;

my $forcerun = "off";
   $forcerun = "on" if defined $opt_F;
   $forcerun = "on" if defined $opt_f;
LOG($outLog, "forcerun is $YW$forcerun$N\n\n");


if (defined $opt_p) {
   LOG($outLog, "#Parallel run with slurm enabled (-p)\n\n");
   if (defined $opt_0) {
      LOG($outLog, date() . "\n" . "$cmd\n");
   }
	my $mem = 16000;
	my $debug;
	$debug = 1 if defined $opt_0;
   my $force_sbatch = $opt_F;
      $force_sbatch = $opt_f if not defined $opt_F;
	footLoop_sbatch_main($sbatch_cmd, "footClust_sbatch", \@neworigFiles, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir, $mem, $debug);
}
else {
   LOG($outLog, "#Single run\n\n");
   for (my $i = 0; $i < @neworigFiles; $i++) {
      my $neworigFile = $neworigFiles[$i];
      my ($neworigFilename) = getFilename($neworigFile, "full");
      my $neworigFileDone = "$outsbatchDir/$neworigFilename.done";
      my $currcmd = $cmd;
      my $indice = $i + 1;
      $currcmd =~ s/FILENAME/$neworigFile/g;
      $currcmd =~ s/FNINDICE/$indice/g;
      if (defined $opt_0) {
         LOG($outLog, date() . "$YW$indice/$totalorigFile$N $LCY$neworigFilename$N\n");
         LOG($outLog, "\t$LGN$currcmd$N\n");
      }
      elsif (not -e $neworigFileDone or defined $opt_F or defined $opt_f) {
         system($currcmd) == 0 or DIELOG($outLog, "Failed to run cmd: $!\n\n$LCY$currcmd$N\n\n");
         system("touch $neworigFileDone") if not -e $neworigFileDone;
      }
      else {
         LOG($outLog, date() . "${LPR}$0$N: $YW$indice/$totalorigFile$N $LCY$neworigFilename$N: using previously made peaks\n");
         LOG($outLog, "Done=$LCY$neworigFileDone$N)\n","NA");
      }
   }
}

LOG($outLog, "\n\n");
LOG($outLog, date() . "${LGN}SUCCESS!!$N\n\n${YW}footClust_sbatch.pl -n $LCY$footPeakFolder$N $LGN$runScript$N\n\n");
LOG($outLog, "\n\n");
