#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_t $opt_R $opt_D $opt_0 $opt_J $opt_F);
getopts("vd:n:G:t:RD:0J:F");

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

my $date = getDate();
my $uuid = getuuid();

my ($dist, $footPeakFolder) = ($opt_d, $opt_n);
my $toggleRstrand = defined $opt_R ? "Yes" : "No";

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N $LGN-g gene$N $CY-n <footPeak's output folder (footPeak's -o)>$N

-D: tightness of read placement in each cluster (default: 200)
-R: toggle to order reads reversely (good foor in vitro REVERSE_ genes)

";

die $usage  unless defined $opt_n and -d $opt_n;

my $outDir = "$footPeakFolder/FOOTCLUST/";
$outDir= getFullpath($outDir);
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
     #    $footPeakFolder/.footClust_sbatch_2/$peakFile\_footClust_logFile.txt

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
my @filesARRAY;
for (my $i = 0; $i < @local_peak_files; $i++) {
	my $local_peak_file = $local_peak_files[$i];
	#run_footClust_sbatch_2(\@cmd);
	my $peakFile = $local_peak_file;
	if (defined $opt_G) {
		if ($local_peak_file !~ /$opt_G/i) {
			#LOG($outLog, date() . "${LPR}footClust_sbatch_2.pl$N: $LGN$i/$totalpeakFiles$N ${LRD}Skipped$N $LCY$local_peak_file$N as it doesn't contain $LGN-G $opt_G$N\n");
			next;
		}
		else {
			LOG($outLog, date() . "${LPR}footClust_sbatch_2.pl$N: $LGN$i/$totalpeakFiles$N ${LGN}Processing$N $LCY$local_peak_file$N as it contain $LGN-G $opt_G$N\n");
		}
	}
	else {
		LOG($outLog, date() . "${LPR}footClust_sbatch_2.pl$N: $LGN$i/$totalpeakFiles$N ${LGN}Processing$N $LCY$local_peak_file$N\n");
	}

	#my @cmd = ("FILENAME", $faFile, $label, $footLoopFolder, $geneIndexFile, $totalpeakFiles);
	push(@filesARRAY, $peakFile);
	#my $cmd = "footClust_sbatch_2.pl";
	#$cmd .= " -n $footPeakFolder -J $max_parallel_run -t $clustThreshold -D $windowDiv -d $dist";
	#$cmd .= " -R" if defined $opt_R;
	#$cmd .= " -F" if defined $opt_F;
	#$cmd .= " -G $opt_G" if defined $opt_G;
	#$cmd .= " " . join(" ", @cmd);
	#LOG($outLog, "Example cmd:\n$LCY$cmd$N\n\n") if $i == 0;
	#LOG($outLog, "$i. footClust_sbatch_2.pl $peakFile\n") if $i > 0;
	#system($cmd) == 0 or die "Failed to run cmd: $!\n";
	#last;
	#last if $i >= 2;
}

my @cmd = ("FILENAME", $faFile, $label, $footLoopFolder, $geneIndexFile, $totalpeakFiles);
my $cmd = "footClust_sbatch_2.pl";
$cmd .= " -n $footPeakFolder";
my $cmd2 = "-J $max_parallel_run -t $clustThreshold -D $windowDiv -d $dist";
$cmd2 .= " -R" if defined $opt_R;
$cmd2 .= " -F" if defined $opt_F;
$cmd2 .= " -G $opt_G" if defined $opt_G;
my $runScript = $cmd2;
$cmd .= " $cmd2 " . join(" ", @cmd);
#LOG($outLog, "Example cmd:\n$LCY$cmd$N\n\n");
my $sbatch_these_cmd = $cmd; #"footClust_sbatch_2.pl -n FILENAME -J ..etc"
my $force_sbatch = 1 if defined $opt_F;
my $outsbatchDir = "$footPeakFolder/.footClust_sbatch/";
system("mkdir -p $outsbatchDir") if not -d $outsbatchDir;
my $mem = 8000;
my $debug;
$debug = 0 if defined $opt_0;
my $force_sbatch_print = defined $force_sbatch ? "force_sbatch is on" : "force_sbatch_is_off";
LOG($outLog, "force sbatch is $YW$force_sbatch_print$N\n\n");
footLoop_sbatch_main($sbatch_these_cmd, "footClust_sbatch", \@filesARRAY, $max_parallel_run, $outLog, $force_sbatch, $outsbatchDir, $mem, $debug);

LOG($outLog, "\n\n");
LOG($outLog, date() . "${LGN}SUCCESS!!$N\n\n${YW}footClust_sbatch.pl -n $LCY$footPeakFolder$N $LGN$runScript$N\n\n");
LOG($outLog, "\n\n");

sub parse_footLoop_logFile {
   my ($logFile, $date, $uuid, $footFolder, $version) = @_;
   #my @line = `cat $logFile`;
   my $paramsFile = "$footFolder/.PARAMS";
 #  my $inputFolder = $defOpts->{n};
	my $geneIndexFile;
   my @parline = `cat $paramsFile`;
   foreach my $parline (@parline) {
      if ($parline =~ /footLoop.pl,geneIndexFile,/) {
         ($geneIndexFile) = $parline =~ /geneIndexFile,(.+)$/;
      }
	}
	return($geneIndexFile);
}

__END__
