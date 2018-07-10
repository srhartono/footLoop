#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Time::HiRes; use Benchmark qw(:all); use Benchmark ':hireswallclock'; use Carp;
use Thread; use Thread::Queue;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_n $opt_i $opt_V $opt_a $opt_b $opt_f $opt_h $opt_H $opt_o $opt_l);
getopts("vn:i:Va:b:f:hHo:l:");
srand(420);
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';

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

# Specify Date, uuid, folder, etc
my $date = getDate();
my $uuid = getuuid();
my ($footFolder) = $opt_n;
my ($scriptName) = getFilename($0);
#my ($FOLDER1) = @ARGV;
my $LABEL = defined $opt_l ? $opt_l : "";
my $usage = "

$YW ------------------------------------------- $N
$YW|$N $scriptName version $version $YW|$N
$YW ------------------------------------------- $N
Script Path: $0

Usage: $YW$scriptName$N $LGN-o$N <outputdir> $LPR-i$N <geneIndexes.bed> $LCY-a$N <bedfile1/folder1> [${LGN}Optional$N: $LGN-b$N <bedfile2/folder2>]

$LGN-o$N: output directory
$LPR-i$N: original footLoop geneIndexes$LCY bed6$N file
$LCY-a$N: footPeak peak bedfile or folder containing peak bedfiles (peaks has to come ${LCY}from PEAKS_GENOME$N)
$LGN-b$N: (Optional): bedfile or folder containing bedfiles to be compared with -a
-l: (optional) label for the output result pdf

Do $YW-h$N for more information about options and file formats
Do $YW-H$N for even longer information

";

my $usage_med = "
This script will group peaks with:
1. same gene
2. same strand
3. same conversion type (e.g. C->T)
Then perform comparison of peak bed files, categorized by:
1. $LGN label$N (e.g. PCB1 or PCB12)
2. $LGN window$N (e.g. 10 or 20)
3. $LGN threshold$N (e.g. 0.35, 0.65)
4. $LGN CpG vs. no CpG$N (e.g. CG vs CH)

${LGN}Options:$N
${YW}1.$N If -a is a folder containing multiple bedfiles and -b not given:
          It'll perform all vs. all comparison between peak bedfiles inside -a folder
${YW}2.$N If -a is a folder containing multiple bedfiles and -b is also a folder with multiple bedfiles:
          It'll perform comparison of all peak bedfiles inside -a folder vs. all peak bedfiles inside -b folder
${YW}3.$N If -a is a$LPR bedfile$N and -b is a folder with multiple bedfiles:
          It'll perform comparison of -a bedfile with all peaks inside -b folder
${YW}4.$N If -a is a$LPR bedfile$N and -b is a also a$LPR bedfile$N:
          It'll perform comparison of -a bedfile with -b bedfile

${LRD}IMPORTANT!!!$N

1. ${LCY}bed files are BED6 from footPeak.pl PEAKS_GENOME folder$N
2. Make sure that the name of bedfile(s) to be compared follows standard footPeak format or it'll bug!
";

my $usage_long = "
   Good example:$LCY PCB5_geneMRFAP1L1_Pos_20_0.55_CH.PEAK.genome.bed$N

   footPeak format:
   - In Vivo data:  <label>_gene<gene>_<strand>_<window>_<threshold>_<type>.PEAK.genome.bed
   - In Vitro data: <label>_<barcode>_gene<gene>_<strand>_<window>_<threshold>_<type>.PEAK.genome.bed


   label=PCB<digits> e.g. PCB5
   barcode=BC<digits>_<alphanumeric description> e.g. BC9
   strand=<Pos|Neg|Unk> e.g. Pos
   window=<integer> e.g. 20
   threshold=<decimal between 0 and 1> e.g. 0.55
   type=<CH|CG|GH|GC> e.g. CH

   Bad examples: - MRFAP1L1_Pos_CH (no label and no window/threshold)
                 - PCB5_geneMRFAP1L1_pos_20_65_CH (strand should be Pos not pos, threshold is not in decimal between 0 to 1)


${YW}Example$N

- Folder1 contains:
1.1. PCB1_geneCALM3_Pos_20_0.65_CH.PEAK.genome.bed$LGN #CALM3 gene from PCB1. Good file$N
1.2. PCB2_geneCALM3_Pos_20_0.65_CH.PEAK.genome.bed$LGN #same as (1) but PCB2 and extension differs$N
1.3. PCB2_geneCALM3_Pos_20_0.35_CH.PEAK.genome.bed$LGN #same as (1) but 0.35 threshold$N and window is 10$N
1.4. PCB1_geneCALM3_Neg_10_0.65_CH.PEAK.genome.bed$LGN #non R-loop forming strand of (1)$N
1.5. PCB5_geneCALM3_Pos_20_0.65_CG.PEAK.genome.bed$LGN #same as (1) but includes CpG$N
1.6. PCB5_geneCALM3_Pos_20_0.65_GC.PEAK.genome.bed$LGN #same as (1) but GC version$N
1.7. PCB5_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed$LGN #FUS gene from PCB5$N
1.8. PCB8_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed$LGN #FUS gene same as (7) but from PCB8$N
1.9. geneCALM3_Pos_20_0.65 #same as (1) but label and conversion type is missing$N

- Folder2 contains:
2.1. PCB3_geneCALM3_Pos_20_0.55_CH.PEAK.genome.bed$LGN #CALM3 pos from PCB3 and 55% threshold$N
2.2. PCB5_geneCALM3_Pos_20_0.55_CH.PEAK.genome.bed$LGN #CALM3 pos from PCB5 and 55% threshold$N
2.3. PCB3_geneFUS_Pos_20_0.55_CH.PEAK.genome.bed$LGN #FUS pos from PCB3 and 55% threshold$N
2.4. PCB3_geneCALM3_Neg_20_0.55_CH.PEAK.genome.bed$LGN #CALM3 neg from PCB3 and 55% threshold$N


${LCY}Example 1:$YW $scriptName$N -a$LGN Folder1$N
-> will group files into these groups and compare them)
- group1: (CALM3 Pos C->T): file 1.1 vs 1.2 vs 1.3 vs 1.5 (4x4 XY plots. Notice that CH and CG (1.5) are compared together)
- group2: (CALM3 Neg C->T): file 1.4 (but only have 1 sample so ignored)
- group3: (CALM3 Pos G->A): file 1.6 (again ignored)
- group4: (FUS Pos C->T)  : file 1.7 and 1.8 
${LRD}1.9 is ignored as sample name doesn't follow standard footPeak format$N

${LCY}Example 2:$YW $scriptName$N -a$LGN Folder1$N -b$LGN Folder2$N:
-> will no longer compare files within each folder like before, but only between relevant files in folder 1 vs. folder2.
- group1: (CALM3 Pos C->T): 1.1 vs. 2.1, 1.1 vs 2.2, 1.2 vs 2.1, 1.2 vs 2.2, etc$LGN but NOT 1.1 vs 1.2 for example$N
- group2: (FUS Pos C->T)  : 1.7 vs 2.3 and 1.8 vs. 2.3$LGN but NOT 1.7 vs 1.8$N
- group3: (CALM3 Neg C->T): 1.4 vs. 2.4
- group4: (CALM3 Pos G->A): 1.6 (is ignored as no other files to compare with)

${LCY}Example 3:$YW $scriptName$N -a$LGN Folder1/PCB1_geneCALM3_Pos_20_0.65_CH.PEAK.genome.bed$N -b$LGN Folder2$N:
-> will only compare file specified in -a vs all relevant files in -b folder2
group1: (CALM3 Pos C->T): 1.1 vs. 2.1 and 1.1 vs 2.2 (1x2 XY plot)

${LCY}Example 4:$YW $scriptName$N -a$LGN Folder1/PCB5_geneCALM3_Pos_20_0.65_GC$N -b$LGN Folder2$N:
-> nothing will be made as it has no files to compare to in folder 2

";

die $usage . "$YW -------- ADDITIONAL INFORMATION -------- $N\n" . $usage_med . "\n\n$YW ----------------------------------------- $N\n\n" if defined $opt_h;
die $usage . "$YW -------- ADDITIONAL INFORMATION -------- $N\n" . $usage_med . "\n\n$YW -------- FILE FORMAT AND EXAMPLES -------- $N" . $usage_long . "$YW ---------------------------------------- $N\n\n" if defined $opt_H;
die $usage . "\n$YW ---------------------------------------- $N\n\n" unless (defined $opt_a and defined $opt_b and -e $opt_a and -e $opt_b) or (defined $opt_a and -e $opt_a);
die $usage . "\n\n$LRD ERROR$N: Please define output dir (-o)\n\n" if not defined $opt_o;
die $usage . "\n\n$LRD ERROR$N: Please input geneIndex (-i)\n\n" if not defined $opt_i;
die $usage . "\n\n$LRD ERROR$N: geneIndex (-i $opt_i) does not exists!\n\n" if not -e $opt_i;
my $outdir = $opt_o;
makedir($outdir) if not -d $outdir;

my ($FOLDER1, $FOLDER2) = ($opt_a, $opt_b);
$FOLDER1 = getFullpath($FOLDER1);
$FOLDER2 = getFullpath($FOLDER2) if defined $opt_b and -e $opt_b;

my $logFile = "$outdir/$scriptName\_logFile.txt";
open (my $outLog, ">", $logFile) or print date() . "\n\n${LRD}ERROR$N: Failed to write to $LCY$logFile$N: $!\n\n" and die;

# Parse gene Index
my $indexFile = defined $opt_i ? $opt_i : "/home/mitochi/pacbio_indexes/0_all_indexes_Alpha301.bed";
LOG($outLog,  "\n" . date() . "${YW}1$N. Parsing geneIndex file -i $LGN$indexFile$N\n");
my $index = parse_indexFile($indexFile);
my $totalgene = (keys %{$index});
LOG($outLog,  date() . "  Parsed $LGN$totalgene$N genes!\n");
# Parse files
my @FOLDER1 = -d $FOLDER1 ? <$FOLDER1/*.PEAK.genome.bed> : $FOLDER1 =~ /.PEAK.genome.bed$/ ? ($FOLDER1) : DIELOG($outLog, date() . "ERROR: -a $FOLDER1 does not end with .PEAK.genome.bed!\n");
my @FOLDER2 = (not defined $opt_b) ? @FOLDER1 : (defined $opt_b and not -e $opt_b) ? DIELOG($outLog, date() . " ERROR There's no gene in $opt_b!\n") : -d $FOLDER2 ? <$FOLDER2/*.PEAK.genome.bed> : $FOLDER2 =~ /.PEAK.genome.bed$/ ? ($FOLDER2) : DIELOG($outLog, date() . "ERROR: -a $FOLDER2 does not end with .PEAK.genome.bed!\n");
my $totFOLDER1 = @FOLDER1;
my $totFOLDER2 = @FOLDER2;
LOG($outLog,  "\n" . date() . "${YW}2$N. Doing all-vs-all peak comparison in -a $LCY$FOLDER1$N ($LGN$totFOLDER1$N .PEAK.local.genome files)\n") if not defined $opt_b;
LOG($outLog,  "\n" . date() . "${YW}2$N. Comparing peaks in -a $LCY$FOLDER1$N ($LGN$totFOLDER1$N .PEAK.local.genome files) vs. -b $LCY$FOLDER2$N ($LGN$totFOLDER2$N .PEAK.local.genome files)\n") if defined $opt_b;

LOG($outLog, "\n\n$YW -------------- Parsing Bed Files from -a $FOLDER1 ----------- $N\n\n");
my %processed;
my ($fileCount, $data, $done) = (0);
my %input;
foreach my $input1 (@FOLDER1[0..@FOLDER1-1]) {
	$fileCount ++;
	my $total = 0;
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	my ($label, $barcode, $treat, $gene, $strand, $window, $thres, $type) = parseNameDetail($fileName1);
	my $conv = $type =~ /^C/ ? "C" : "G";
	DIELOG($outLog, date() . "$LRD  ERROR$N: gene $gene doesn't exists in indexFile $LCY$indexFile$N\n") if not defined $index->{$gene};
	($data, $total) = parse_bedFile($input1, $data);
	LOG($outLog, date() . "$YW$fileCount$N. Parsed $LGN$total$N peaks from label=$label,bc=$barcode,treat=$treat,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type,conv=$conv,file=$LCY$fileName1$N\n");
#	LOG($outLog, date() . "\t- Parsed $LGN$total$N peaks\n");
	$input{"gene$gene\_strand$strand\_conv$conv"}{1}{$input1} = 1;
	$input{"gene$gene\_strand$strand\_conv$conv"}{2}{$input1} = 1 if not defined $opt_b;
}
if (defined $opt_b and -e $opt_b) {
	my $fileCount2 = 0;
	LOG($outLog, "\n\n$YW -------------- Parsing Bed Files from -a $FOLDER2 ----------- $N\n\n") if defined $opt_b;
	foreach my $input2 (@FOLDER2[0..@FOLDER2-1]) {
		$fileCount ++;
		$fileCount2 ++;
		my $total = 0;
		my ($folder2, $fileName2) = getFilename($input2, "folderfull");
		my ($label, $barcode, $treat, $gene, $strand, $window, $thres, $type) = parseNameDetail($fileName2);
		my $conv = $type =~ /^C/ ? "C" : "G";
		DIELOG($outLog, date() . "$LRD  ERROR$N: gene $gene doesn't exists in indexFile $LCY$indexFile$N\n") if not defined $index->{$gene};
		($data, $total) = parse_bedFile($input2, $data);
		LOG($outLog, date() . "$YW$fileCount2$N. Parsed $LGN$total$N peaks from label=$label,bc=$barcode,treat=$treat,gene=$gene,strand=$strand,window=$window,thres=$thres,type=$type,conv=$conv,file=$LCY$fileName2$N\n");
	#	LOG($outLog, date() . "$YW$fileCount$N. Parsed $LGN$total$N from $LCY$input2$N\n");
		$input{"gene$gene\_strand$strand\_conv$conv"}{2}{$input2} = 1;
	}
}

foreach my $group (sort keys %input) {
	print "\n> Group $YW$group$N\n";
	my $file1count = 0;
	foreach my $input1 (sort keys %{$input{$group}{1}}) {
		my ($fileName1) = getFilename($input1, "full");
		$file1count ++;
		my $file2count = 0;
		foreach my $input2 (sort keys %{$input{$group}{2}}) {
			my ($fileName2) = getFilename($input2, "full");
			next if $input1 eq $input2;
			$file2count ++;
			print "$file1count vs $file2count: $LGN$fileName1$N vs $LCY$fileName2$N\n" if not defined $processed{$group}{"$input2\_$input1"} and not defined $processed{$group}{"$input1\_$input2"};
			$processed{$group}{"$input2\_$input1"} = 0;
			$processed{$group}{"$input1\_$input2"} = 0;
		}
	}
}
LOG($outLog, "\n\n$YW -------------- Calculating Distances ----------- $N\n\n");

$fileCount = 0;
my $geneCount = 0;
open (my $out1, ">", "$outdir/RESULTS.tsv") or DIELOG($outLog, "\n\nFailed to write to $LCY$FOLDER1/RESULTS.tsv$N: $!\n\n");
#print $out1 "shuftype\trestype\tsample1\ttype1\tsample2\ttype2\tstrand\tgene\tpeak1\tpeak2\tdiff\txpos\typos\tbeg1\tend1\tbeg2\tend2\tmisc\tpeakname\n";
print $out1 "shuftype\trestype1\trestype2\tgroup\tpeak1\tpeak2\tdiff\txpos\typos\tchr1\tbeg1\tend1\tchr2\tbeg2\tend2\tmisc\tpeakname\n";
foreach my $group (sort keys %input) {
	my ($genez, $strandz, $convz) = $group =~ /^gene(.+)_strand(.+)_conv(.+)$/;
	DIELOG($outLog, date() . " ERROR undefined genez=$genez/strandz=$strandz/convz=$convz from group=$group\n\n") if not defined $genez or not defined $strandz or not defined $convz;
#	next unless $genez eq "CALM3";
	LOG($outLog, "\n" . date() . "$YW$geneCount$N. Processing group: $LCY$group$N\n");
		foreach my $input1 (sort keys %{$input{$group}{1}}) {
		$fileCount ++;
	
		LOG($outLog, date() . " ...$YW$geneCount$N.$LGN$fileCount$N. Processing input1: $LCY$input1$N\n");
		my ($folder1, $fileName1) = getFilename($input1, "folderfull");
		my ($label1, $barcode1, $treat1, $gene1, $strand1, $window1, $thres1, $type1) = parseNameDetail($fileName1);
		my $conv1 = $type1 =~ /^C/ ? "C" : "G";
		my $labelz1 = "$label1\_$window1\_$thres1\_$type1" if $barcode1 eq "";
		   $labelz1 = "$label1\_$barcode1\_$treat1\_$window1\_$thres1\_$type1" if $barcode1 ne "";
		my $fileCount2 = 0;
		foreach my $input2 (sort keys %{$input{$group}{2}}) {
			#next unless $input1 =~ /PCB1_/ and $input2 =~ /PCB5_/ and $genez eq "CALM3";
			next if defined $done->{$group}{"$input1\_$input2"} or defined $done->{$group}{"$input2\_$input1"};
			$fileCount2 ++; #last if $fileCount2 > 2;
			my $STRAND = $index->{$genez}{strand};
			$STRAND = "Pos" if $STRAND eq "+";
			$STRAND = "Neg" if $STRAND eq "-";
#			LOG($outLog, date() . " .....$YW$geneCount$N.$LGN$fileCount$N.$LCY$fileCount2$N. vs. input2: $LCY$input2$N\n");
			my ($folder2, $fileName2) = getFilename($input2, "folderfull");
			my ($label2, $barcode2, $treat2, $gene2, $strand2, $window2, $thres2, $type2) = parseNameDetail($fileName2);
		
			my $conv2 = $type2 =~ /^C/ ? "C" : "G";
			my $labelz2 = "$label2\_$window2\_$thres2\_$type2" if $barcode2 eq "";
			   $labelz2 = "$label2\_$barcode2\_$treat2\_$window2\_$thres2\_$type2" if $barcode2 ne "";

			$done->{$group}{"$input1\_$input2"} = 1;
			$done->{$group}{"$input2\_$input1"} = 1;
			$processed{$group}{"$input1\_$input2"} = 1;
			$processed{$group}{"$input2\_$input1"} = 1;
#			next if ($gene1 ne $gene2 or $strand1 ne $strand2 or $conv or $window1 ne $window2 or $thres1 ne $thres2 or $type1 ne $type2);
			#my ($label1, $gene1, $strand, $window, $thres, $type) = ($label1, $gene1, $strand1, $window1, $thres1, $type1);
	
			my $restype1 = "OTHERS";
			my $restype2 = "OTHERS";
			if ($label1 !~ /PCB(11|13|14|15|17)/) {
				$restype1 = (($strand1 eq "Pos" and $type1 eq "CH" and $STRAND eq "Pos") or ($strand1 eq "Neg" and $type1 eq "GH" and $STRAND eq "Neg")) ? "PEAK" :
							   (($strand1 eq "Pos" and $type1 eq "CH" and $STRAND eq "Neg") or ($strand1 eq "Neg" and $type1 eq "GH" and $STRAND eq "Pos")) ? "PEAKNEG": $restype1;
			}
			else {
				$restype1 = (($strand1 eq "Pos" and $type1 eq "CG" and $STRAND eq "Pos") or ($strand1 eq "Neg" and $type1 eq "GC" and $STRAND eq "Neg")) ? "PEAK" :
							   (($strand1 eq "Pos" and $type1 eq "CG" and $STRAND eq "Neg") or ($strand1 eq "Neg" and $type1 eq "GC" and $STRAND eq "Pos")) ? "PEAKNEG": $restype1;
			}
			if ($label2 !~ /PCB(11|13|14|15|17)/) {
				$restype2 = (($strand2 eq "Pos" and $type2 eq "CH" and $STRAND eq "Pos") or ($strand2 eq "Neg" and $type2 eq "GH" and $STRAND eq "Neg")) ? "PEAK" :
							   (($strand2 eq "Pos" and $type2 eq "CH" and $STRAND eq "Neg") or ($strand2 eq "Neg" and $type2 eq "GH" and $STRAND eq "Pos")) ? "PEAKNEG": $restype1;
			}
			else {
				$restype2 = (($strand2 eq "Pos" and $type2 eq "CG" and $STRAND eq "Pos") or ($strand2 eq "Neg" and $type2 eq "GC" and $STRAND eq "Neg")) ? "PEAK" :
							   (($strand2 eq "Pos" and $type2 eq "CG" and $STRAND eq "Neg") or ($strand2 eq "Neg" and $type2 eq "GC" and $STRAND eq "Pos")) ? "PEAKNEG": $restype1;
			}
			my ($DRIPtype1, $DRIPtype2) = ("DRIP","DRIP");
			$DRIPtype1 = "noDRIP" if ($input1 =~ /PCB(2|8|9|10|11|13|14|15|16|17)_/);
			$DRIPtype2 = "noDRIP" if ($input2 =~ /PCB(2|8|9|10|11|13|14|15|16|17)_/);
			$restype1 = "$labelz1\_$restype1\_$DRIPtype1";
			$restype2 = "$labelz2\_$restype2\_$DRIPtype2";
			
			LOG($outLog, date() . "   -->$YW$geneCount$N.$LGN$fileCount$N.$LPR$fileCount2$N. Processing $LCY$input1$N restype1=$LGN$restype1$N restype2=$LPR$restype2$N\n");

			my @orig1 = @{$data->{$input1}};
			my @orig2 = @{$data->{$input2}};
			my $BEG = $index->{$genez}{beg};
			my $END = $index->{$genez}{end};
			my $added1 = 0;
			my $added2 = 0;
			LOG($outLog, date() . "\t\ttotal:peak1=" . scalar(@orig1) . ",peak2=" . scalar(@orig2) . ",");
			for (my $i = 0; $i < @orig1; $i++) {
				my ($chr, $beg, $end) = split(",", $orig1[$i]);
				$beg -= $BEG;
				$end -= $BEG;
				$orig1[$i] = "$chr,$beg,$end";
			}
			for (my $i = 0; $i < @orig2; $i++) {
				my ($chr, $beg, $end) = split(",", $orig2[$i]);
				$beg -= $BEG;
				$end -= $BEG;
				$orig2[$i] = "$chr,$beg,$end";
			}
			for (my $i = @orig1; $i < @orig2; $i++) {
				my $rand = $orig1[int(rand(@orig1))];
				$orig1[$i] = $rand;
				$added1 ++;
			}
			for (my $i = @orig2; $i < @orig1; $i++) {
				my $rand = $orig2[int(rand(@orig2))];
				$orig2[$i] = $rand;
				$added2 ++;
			}
			LOG($outLog, " added $LGN$added1$N to peak1 and $LGN$added2$N to peak2\n");

			my %orig;
			my %shuf1; my %shuf2;
			my $LEN = $END - $BEG;
			my %orig3; my %shuf3;
			my %orig4; my %shuf4;
			for (my $i = 0; $i < @orig1; $i++) {
				my ($chr, $beg, $end) = split(",", $orig1[$i]);
				$orig3{$i}{chr} = $chr;
				$orig3{$i}{beg} = $beg;
				$orig3{$i}{end} = $end;
				$orig3{$i}{begindex} = int($beg/100);
				$orig3{$i}{endindex} = int($end/100);
			}
			for (my $i = 0; $i < @orig2; $i++) {
				my ($chr, $beg, $end) = split(",", $orig2[$i]);
				$orig4{$i}{chr} = $chr;
				$orig4{$i}{beg} = $beg;
				$orig4{$i}{end} = $end;
				$orig4{$i}{begindex} = int($beg/100);
				$orig4{$i}{endindex} = int($end/100);
			}
			my @shuf1; my @shuf2;
			for (my $i = 0; $i < @orig1; $i ++) {
				my ($chr, $beg, $end) = split(",", $orig1[$i]);
				my $len = $end - $beg;
				my $lendiff = $LEN - $len;
				my $begrand = int(rand($lendiff));# + $LEN;
				my $endrand = $begrand + $len;# + $LEN;
#				$beg += $LEN;
#				$end += $LEN;
				$shuf1[$i] = "$chr,$begrand,$endrand";
				$shuf3{$i}{chr} = $chr;
				$shuf3{$i}{beg} = $begrand;
				$shuf3{$i}{end} = $endrand;
				$shuf3{$i}{begindex} = int($begrand/100);
				$shuf3{$i}{endindex} = int($endrand/100);
#				print "BEG=$BEG, END=$END, LEN=$LEN, beg=$beg, end=$end, len=$len, shuf1 = begrand=int(rand($lendiff)), endrand = $begrand + $len = $begrand,$endrand\n" if $i < 20;
			}
			for (my $i = 0; $i < @orig2; $i ++) {
				my ($chr, $beg, $end) = split(",", $orig2[$i]);
				my $len = $end - $beg;
				my $lendiff = $LEN - $len;
				my $begrand = int(rand($lendiff));# + $LEN;
				my $endrand = $begrand + $len;# + $LEN;
				$shuf2[$i] = "$chr,$begrand,$endrand";
				$shuf4{$i}{chr} = $chr;
				$shuf4{$i}{beg} = $begrand;
				$shuf4{$i}{end} = $endrand;
				$shuf4{$i}{begindex} = int($begrand/100);
				$shuf4{$i}{endindex} = int($endrand/100);
#				print "BEG=$BEG, END=$END, LEN=$LEN, beg=$beg, end=$end, len=$len, shuf2 = begrand=int(rand($lendiff)), endrand = $begrand + $len = $begrand,$endrand\n" if $i < 20;
			}
			my (@orig3, @orig4, @shuf3, @shuf4);
			#@orig3 = sort {$orig3[$a][2] <=> $orig3[$b][2] || $orig3[$a][3] <=> $orig3[$b][3]} @orig3;
			#@orig4 = sort {$orig4[$a][2] <=> $orig4[$b][2] || $orig4[$a][3] <=> $orig4[$b][3]} @orig4;
			#@shuf3 = sort {$shuf3[$a][2] <=> $shuf3[$b][2] || $shuf3[$a][3] <=> $shuf3[$b][3]} @shuf3;
			#@shuf4 = sort {$shuf4[$a][2] <=> $shuf4[$b][2] || $shuf4[$a][3] <=> $shuf4[$b][3]} @shuf4;
			foreach my $num (sort {$orig3{$a}{begindex} <=> $orig3{$b}{begindex} || $orig3{$a}{endindex} <=> $orig3{$b}{endindex}} keys %orig3) {	
				my $chr = $orig3{$num}{chr};
				my $beg = $orig3{$num}{beg};
				my $end = $orig3{$num}{end};
				push(@orig3, "$chr,$beg,$end");
			}
			foreach my $num (sort {$orig4{$a}{begindex} <=> $orig4{$b}{begindex} || $orig4{$a}{endindex} <=> $orig4{$b}{endindex}} keys %orig4) {	
				my $chr = $orig4{$num}{chr};
				my $beg = $orig4{$num}{beg};
				my $end = $orig4{$num}{end};
				push(@orig4, "$chr,$beg,$end");
			}
			foreach my $num (sort {$shuf3{$a}{begindex} <=> $shuf3{$b}{begindex} || $shuf3{$a}{endindex} <=> $shuf3{$b}{endindex}} keys %shuf3) {	
				my $chr = $shuf3{$num}{chr};
				my $beg = $shuf3{$num}{beg};
				my $end = $shuf3{$num}{end};
				push(@shuf3, "$chr,$beg,$end");
			}
			foreach my $num (sort {$shuf4{$a}{begindex} <=> $shuf4{$b}{begindex} || $shuf4{$a}{endindex} <=> $shuf4{$b}{endindex}} keys %shuf4) {	
				my $chr = $shuf4{$num}{chr};
				my $beg = $shuf4{$num}{beg};
				my $end = $shuf4{$num}{end};
				push(@shuf4, "$chr,$beg,$end");
			}	
			for (my $i = 0; $i < @orig3; $i++) {
				my ($chr1, $beg1, $end1) = split(",", $orig3[$i]);
				my ($chr2, $beg2, $end2) = split(",", $orig4[$i]);
				my ($xpos, $ypos) = ($end2-$beg1+$LEN, $end1-$beg2+$LEN);
				my $coors = "$chr1\t$beg1\t$end1\t$chr2\t$beg2\t$end2";
#				print $out1 "ORIG\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak1\n";# if $peaktemp == 2;
#				print $out1 "ORIG\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak2\n";# if $peaktemp == 2;
			}
			for (my $i = 0; $i < @orig3; $i++) {
				my ($chr1, $beg1, $end1) = split(",", $orig3[$i]);
				my ($chr2, $beg2, $end2) = split(",", $shuf4[$i]);
				my ($xpos, $ypos) = ($end2-$beg1+$LEN, $end1-$beg2+$LEN);
				my $coors = "$chr1\t$beg1\t$end1\t$chr2\t$beg2\t$end2";
#				print $out1 "SHUF1\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak1\n";# if $peaktemp == 2;
#				print $out1 "SHUF1\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak2\n";# if $peaktemp == 2;
			}
			for (my $i = 0; $i < @orig3; $i++) {
				my ($chr1, $beg1, $end1) = split(",", $shuf3[$i]);
				my ($chr2, $beg2, $end2) = split(",", $orig4[$i]);
				my ($xpos, $ypos) = ($end2-$beg1+$LEN, $end1-$beg2+$LEN);
				my $coors = "$chr1\t$beg1\t$end1\t$chr2\t$beg2\t$end2";
#				print $out1 "SHUF2\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak1\n";# if $peaktemp == 2;
#				print $out1 "SHUF2\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\tpeak1\tpeak2\tdiff\t$xpos\t$ypos\t$coors\tmisc\tpeak2\n";# if $peaktemp == 2;
			}
#			next;
			#orig1 and orig2
			my @shuf = qw(ORIG SHUF1 SHUF2);
			#my @peaks1 = (\@orig1, \@orig1, \@shuf1);
			#my @peaks2 = (\@orig2, \@shuf2, \@orig2);
			for (my $h = 0; $h < 3; $h++) {
				my @peak1 = $h == 0 ? @orig1 : $h == 1 ? @orig1 : @shuf1; #@{$peaks1[$h]};# > @{$peaks2[$h]} ? @{$peaks1[$h]} : @{$peaks2[$h]};
				my @peak2 = $h == 0 ? @orig2 : $h == 1 ? @shuf2 : @orig2; #@{$peaks1[$h]};# > @{$peaks2[$h]} ? @{$peaks1[$h]} : @{$peaks2[$h]}; my @peak2 = @{$peaks2[$h]};# <= @{$peaks2[$h]} ? @{$peaks1[$h]} : @{$peaks2[$h]};
				my $peaktemp = 1;#@{$peaks1[$h]} > @{$peaks2[$h]} ? 1 : 2;
				my %res;
				for (my $i = 0; $i < @peak1; $i ++) {
					for (my $j = 0; $j < @peak2; $j ++) {
						my ($chr1, $beg1, $end1) = split(",", $peak1[$i]);
						my ($chr2, $beg2, $end2) = split(",", $peak2[$j]);
						my ($begdiff, $enddiff) = ($beg2-$beg1, $end2-$end1);
						my ($len1, $len2)       = ($end1-$beg1, $end2-$beg2);
						my ($xpos, $ypos)       = ($end2-$beg1, $end1-$beg2);
						my $diff = ($begdiff**2 + $enddiff**2) / ($len1+$len2);
						($diff) = $diff =~ /^(.+\.[0]*[1-9]\d)/ if $diff =~ /^.+\.[0]*[1-9]\d/;
						my $peakpos = "$i,$j";
						$res{$peakpos}{diff} = $diff;
						$res{$peakpos}{peak1} = $i;
						$res{$peakpos}{peak2} = $j;
						$res{$peakpos}{xpos} = $xpos + $LEN;
						$res{$peakpos}{ypos} = $ypos + $LEN;
						$res{$peakpos}{coor} = "$chr1\t$beg1\t$end1\t$chr2\t$beg2\t$end2";
#						$res{$peakpos}{coor} = "$beg2\t$end2\t$beg1\t$end1";# if $peaktemp == 2;
						$res{$peakpos}{misc} = "len1=$len1,len2=$len2,begdiff=$begdiff,enddiff=$enddiff";
						$res{$peakpos}{misc} = "len1=$len2,len2=$len1,begdiff=$begdiff,enddiff=$enddiff" if $peaktemp == 2;
					}			
				}
	
				# peak1 > paek2
				my (%best, %used);
				my $fraction = 2;#@peak1 == @peak2 ? 2 : int(@peak1 / @peak2) + 1;
				my $max = 9999999999;#@peak1 - @peak2;
				my $usedtot = 0;
				foreach my $peakpos (sort {$res{$a}{diff} <=> $res{$b}{diff} || $res{$a}{peak1} <=> $res{$b}{peak1} || $res{$a}{peak2} <=> $res{$b}{peak2}} keys %res) {
					my ($peak1, $peak2) = split(",", $peakpos);
					next if defined $used{"1"}{$peak1};
#					next if defined $used{"1b"}{$peak2};
					if (defined $used{"1b"}{$peak2}) {
						next if $used{"1b"}{$peak2} >= $fraction;
#						next if $usedtot >= $max;
#						$usedtot ++;
					}
					$usedtot ++;
					$best{"1"}{$peak1}{peak} = $peak2;
					$best{"1"}{$peak1}{diff} = $res{$peakpos}{diff};
					$best{"1"}{$peak1}{xpos} = $res{$peakpos}{xpos};
					$best{"1"}{$peak1}{ypos} = $res{$peakpos}{ypos};
					$best{"1"}{$peak1}{coor} = $res{$peakpos}{coor};
					$best{"1"}{$peak1}{misc} = $res{$peakpos}{misc};
					$used{"1"}{$peak1} = 1;
					$used{"1b"}{$peak2} ++;
				}
				my $used1tot = (keys %{$used{"1"}}); my $peak1tot = @peak1;
				DIELOG($outLog, date() . "peaktemp=$peaktemp. Used1 $LGN$used1tot$N is not the same as number of peak1 $LGN$peak1tot$N!\n\n") if $used1tot != $peak1tot;

				$usedtot = 0;
				$fraction = 2;#@peak1 == @peak2 ? 2 : int(@peak2 / @peak1) + 1;
				foreach my $peakpos (sort {$res{$a}{diff} <=> $res{$b}{diff} || $res{$a}{peak1} <=> $res{$b}{peak1} || $res{$a}{peak2} <=> $res{$b}{peak2}} keys %res) {
					my ($peak1, $peak2) = split(",", $peakpos);
					next if defined $used{"2"}{$peak2};
#					next if defined $used{"2b"}{$peak1};
					if (defined $used{"2b"}{$peak1}) {
						next if $used{"2b"}{$peak1} >= $fraction;
#						next if $usedtot >= $max;
#						$usedtot ++;
					}
					$usedtot ++;
					$best{"2"}{$peak2}{peak} = $peak1;
					$best{"2"}{$peak2}{diff} = $res{$peakpos}{diff};
					$best{"2"}{$peak2}{xpos} = $res{$peakpos}{xpos};
					$best{"2"}{$peak2}{ypos} = $res{$peakpos}{ypos};
					$best{"2"}{$peak2}{coor} = $res{$peakpos}{coor};
					$best{"2"}{$peak2}{misc} = $res{$peakpos}{misc};
					$used{"2"}{$peak2} = 1;
					$used{"2b"}{$peak1} ++;
				}
	
				my $used2tot = (keys %{$used{"2"}}); my $peak2tot = @peak2;
				DIELOG($outLog, date() . "peaktemp=$peaktemp. Used2 $LGN$used2tot$N is not the same as number of peak2 $LGN$peak2tot$N!\n\n") if $used2tot != $peak2tot;
				foreach my $peak1 (sort {$a <=> $b || $best{"1"}{$a}{peak} <=> $best{"1"}{$b}{peak}} keys %{$best{"1"}}) {
					my $peak2 = $best{"1"}{$peak1}{peak};
					my $xpos = $best{"1"}{$peak1}{xpos};
					my $ypos = $best{"1"}{$peak1}{ypos};
					my $coor = $best{"1"}{$peak1}{coor};
					my $misc = $best{"1"}{$peak1}{misc};
					my $diff = $best{"1"}{$peak1}{diff};
#					print $out1 "$shuf[$h]\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\t$peak1\t$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak1\n" if $peaktemp == 1;
#					print $out1 "$shuf[$h]\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\t$peak1\t$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak2\n";# if $peaktemp == 2;
					print $out1 "$shuf[$h]\t$restype1\t$restype2\t$group\tpeak1.$peak1\tpeak2.$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak1\n";# if $peaktemp == 2;
				}
				foreach my $peak2 (sort {$a <=> $b || $best{"2"}{$a}{peak} <=> $best{"2"}{$b}{peak}} keys %{$best{"2"}}) {
					my $peak1 = $best{"2"}{$peak2}{peak};
					my $xpos = $best{"2"}{$peak2}{xpos};
					my $ypos = $best{"2"}{$peak2}{ypos};
					my $coor = $best{"2"}{$peak2}{coor};
					my $misc = $best{"2"}{$peak2}{misc};
					my $diff = $best{"2"}{$peak2}{diff};
#					print $out1 "$shuf[$h]\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\t$peak1\t$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak2\n" if $peaktemp == 1;
#					print $out1 "$shuf[$h]\t$restype\t$label1\t$DRIPtype1\t$label2\t$DRIPtype2\t$strand\t$gene\t$peak1\t$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak2\n";# if $peaktemp == 2;
					print $out1 "$shuf[$h]\t$restype1\t$restype2\t$group\tpeak1.$peak1\tpeak2.$peak2\t$diff\t$xpos\t$ypos\t$coor\t$misc\tpeak2\n";# if $peaktemp == 2;
				}
			}
	
		}
	}
}

my $Rscript = Rscript($outdir, "RESULTS.tsv");
open (my $outR, ">", "$outdir/RESULTS.R") or DIELOG($outLog, "\n\n" . date() . " Failed to write to $outdir/RESULTS.R: $!\n");
print $outR $Rscript;
close $outR;


LOG($outLog, "\n\n$YW -------------- Checking if All Files were Processed  ----------- $N\n\n");

foreach my $input12 (sort keys %processed) {
#	print "Not done: $input12\n" if $processed{$group}{$input12} == 0;
}

my $RscriptLog = `run_Rscript.pl $outdir/RESULTS.R`;
LOG($outLog, "\n\n$YW ----------- Ran R script ------------- $N\n\n$RscriptLog\n\n");
=comment
		my $difftot = 0;
		for (my $i = 0; $i < @orig1; $i ++) {
			my $orig1 = $i;
			if (not defined $best{"1"}{$i}) {
				$undef1 ++;
				#print "$LGN$orig1$N\t-1\t0\n";
				next;
			}
			my $orig2 = $best{"1"}{$peak1}{peak};
			my $diff = $best{"1"}{$peak1}{diff};
			my $coor = $best{"1"}{$peak1}{coor};
			my $xpos = $best{"1"}{$peak1}{xpos};
			my $ypos = $best{"1"}{$peak1}{ypos};
			#if (defined $used1{$peak1} or defined $used2{$orig2}) {
			#	die "Already defined peak1=$used1{$peak1} or orig2=$used2{$orig2}\n";
			#}
			$used1{$peak1} = 1;
			$used2{$orig2} = 1;
			print $out1 "$label{$input1}\t$type1\t$label{$input2}\t$type2\t$gene\tpeak$peak1\tpeak$orig2\t$diff\t$xpos\t$ypos\t$coor\n";
		}
		print "undef = $undef1\n";
		next;
		print "2:\n";
		for (my $i = 0; $i < @peak1; $i ++) {
			my $peak1 = $i;
			if (not defined $best{"1"}{$i}) {
				#$undef1 ++;
				print "$LGN$peak1$N\t-1\t0\n";
				next;
			}
		}
		print "3:\n";
		for (my $i = 0; $i < @peak2; $i ++) {
			my $peak2 = $i;
			next if defined $used2{$peak2};
			if (not defined $best{2}{$i}) {
				$undef2 ++;
				print "-1\t$LGN$peak2$N\t0\n";
				next;
			}
			my $peak1 = $best{2}{$peak2}{peak};
			my $diff = $best{2}{$peak2}{diff};
			if (defined $used1{$peak1} or defined $used2{$peak2}) {
				die "Already defined peak2=$used1{$peak2} or peak2=$used2{$peak2}\n";
			}
			$used1{$peak1} = 1;
			$used2{$peak2} = 1;
			print "$peak1\t$LGN$peak2$N\t$diff\n";
		}
		print "undef1 = $undef1 / " . scalar(@peak1) . ", undef2 = $undef2 / " . scalar(@peak2) . "\n";
#		die;
	}
}
=cut



sub parseNameDetail {
	my ($fileName1) = @_;
	my ($label, $gene, $strand, $window, $thres, $type) = parseName($fileName1);
	my ($label1, $label2, $label3) = ("", "", ""); my $barcode = ""; my $treat = "";
	if ($label =~ /_BC\d+/) {
		($label3) = $label =~ /BC\d+_(.+)$/ if $label =~ /BC\d+_.+$/;
		($label1, $label2) = $label =~ /^(PCB\d+)_(BC\d+)/;
		DIELOG($outLog, "\n\nFailed to parse label1/2/3 from label=$LCY$label$N which contain _BC\\d+\n\n") if not defined $label1 or not defined $label2 or not defined $label3;
		$label = $label1;
		$barcode = $label2; $treat = $label3;
	}
	my $labelz = ($label2 eq "" and $label3 eq "") ? "l=$label" : "l=$label;l2=$label2;l3=$label3";
	print "$LGN$fileCount$N. $LCY$fileName1$N: $labelz, g=$gene, s=$strand, w=$window, t=$thres, t=$type\n" if defined $opt_V;
	return($label, $barcode, $treat, $gene, $strand, $window, $thres, $type, $label1, $label2, $label3);
}

sub parse_bedFile {
	my ($input1, $data, $outLog) = @_;
	my $total = 0;
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		my ($chr, $beg, $end, $name, $zero, $strand, $peakFile) = split("\t", $line);
		my ($gene, $read) = $name =~ /^(.+)\.(.+)$/;
		DIELOG($outLog, "\n\nFailed to get gene and read from column 1 in file = $LCY$input1$N\n\n") if not defined $gene or not defined $read;
		my $mid = int(($end - $beg)/2+0.5);
		push(@{$data->{$input1}}, "$chr,$beg,$end");
		#$data->{$input1}{$mid}{$beg}{$end} ++;
		$total ++;
	}
	close $in1;
	return ($data, $total);
}

sub Rscript {
	my ($outdir, $result) = @_;
	$LABEL .= "_" if $LABEL ne "";
	my $Rscript = "
.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\", \"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\", .libPaths()) )
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)

setwd(\"$outdir\")
dm = read.table(\"$result\",header=T,sep=\"\\t\")
id = as.data.frame(plyr::count(dm,vars=c(\"restype1\",\"restype2\")))

for (i in 1:dim(id)[1]) {
print(paste(\"Doing\",i,\": \",id\$restype1[i],\" vs. \",id\$restype2[i],sep=\"\"))
df = dm[dm\$restype1 == id\$restype1[i] & dm\$restype2 == id\$restype2[i],]
origA = df[df\$shuftype == \"ORIG\" & df\$peakname == \"peak1\",]
cor(origA\$xpos,origA\$ypos)
origB = df[df\$shuftype == \"ORIG\" & df\$peakname == \"peak2\",]
shuf1A = df[df\$shuftype == \"SHUF1\" & df\$peakname == \"peak1\",]
shuf1B = df[df\$shuftype == \"SHUF1\" & df\$peakname == \"peak2\",]
shuf2A = df[df\$shuftype == \"SHUF2\" & df\$peakname == \"peak1\",]
shuf2B = df[df\$shuftype == \"SHUF2\" & df\$peakname == \"peak2\",]
origA = origA[order(as.integer(origA\$beg1/100),as.integer(origA\$end1/100),as.integer((origA\$end1-origA\$beg1)/100)),];origA\$y = seq(1,dim(origA)[1])
origB = origB[order(as.integer(origB\$beg2/100),as.integer(origB\$end2/100),as.integer((origB\$end2-origB\$beg2)/100)),];origB\$y = seq(1,dim(origB)[1])
shuf1A = shuf1A[order(as.integer(shuf1A\$beg1/100),as.integer(shuf1A\$end1/100),as.integer((shuf1A\$end1-shuf1A\$beg1)/100)),];shuf1A\$y = seq(1,dim(shuf1A)[1])
shuf1B = shuf1B[order(as.integer(shuf1B\$beg2/100),as.integer(shuf1B\$end2/100),as.integer((shuf1B\$end2-shuf1B\$beg2)/100)),];shuf1B\$y = seq(1,dim(shuf1B)[1])
shuf2A = shuf2A[order(as.integer(shuf2A\$beg1/100),as.integer(shuf2A\$end1/100),as.integer((shuf2A\$end1-shuf2A\$beg1)/100)),];shuf2A\$y = seq(1,dim(shuf2A)[1])
shuf2B = shuf2B[order(as.integer(shuf2B\$beg2/100),as.integer(shuf2B\$end2/100),as.integer((shuf2B\$end2-shuf2B\$beg2)/100)),];shuf2B\$y = seq(1,dim(shuf2B)[1])
origA.cor = cor(origA\$ypos,origA\$xpos,method=\"pearson\")
origB.cor = cor(origB\$ypos,origB\$xpos,method=\"pearson\")
shuf1A.cor = cor(shuf1A\$ypos,shuf1A\$xpos,method=\"pearson\")
shuf1B.cor = cor(shuf1B\$ypos,shuf1B\$xpos,method=\"pearson\")
shuf2A.cor = cor(shuf2A\$ypos,shuf2A\$xpos,method=\"pearson\")
shuf2B.cor = cor(shuf2B\$ypos,shuf2B\$xpos,method=\"pearson\")

outpdf = paste($LABEL,id\$restype1[i],id\$restype2[i],sep=\"\")
pdf(outpdf,width=42,height=7);
p1 = ggplot(origA,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(origA\$restype1),hjust=0)
p2 = ggplot(origA,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(origA\$restype2),hjust=0)
p3 = ggplot(origA,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=origA.cor) + theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(origB,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(origB\$restype2),hjust=0)
p5 = ggplot(origB,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(origB\$restype1),hjust=0)
p6 = ggplot(origB,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=origB.cor)+ theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
p1 = ggplot(shuf1A,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(shuf1A\$restype1),hjust=0)
p2 = ggplot(shuf1A,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=paste(unique(shuf1A\$restype2),\"shuf\"),hjust=0)
p3 = ggplot(shuf1A,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=shuf1A.cor) + theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(shuf1B,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=paste(unique(shuf1B\$restype2),\"shuf\"),hjust=0)
p5 = ggplot(shuf1B,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(shuf1B\$restype1),hjust=0)
p6 = ggplot(shuf1B,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=shuf1B.cor) +theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
p1 = ggplot(shuf2A,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=paste(unique(shuf2A\$restype1),\"shuf\"),hjust=0)
p2 = ggplot(shuf2A,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(shuf2A\$restype2),hjust=0)
p3 = ggplot(shuf2A,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=shuf2A.cor) +theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(shuf2B,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=unique(shuf2B\$restype2),hjust=0)
p5 = ggplot(shuf2B,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color=\"black\") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom=\"text\",x=0,y=0,label=paste(unique(shuf2B\$restype1),\"shuf\"),hjust=0)
p6 = ggplot(shuf2B,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom=\"segment\",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method=\"lm\") + annotate(geom=\"text\",x=2000,y=2000,label=shuf2B.cor) +theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
dev.off()


}
";

	return ($Rscript);
}




__END__


my %data;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

__END__
		#orig2 and orig1
		$fraction = int(@orig1 / @orig2) + 1;
		foreach my $peakpos (sort {$orig{$a}{diff} <=> $orig{$b}{diff} || $orig{$a}{peak1} <=> $orig{$b}{peak1} || $orig{$a}{peak2} <=> $orig{$b}{peak2}} keys %orig) {
			my ($orig1, $orig2) = split(",", $peakpos);
			my $diff = $orig{$peakpos}{diff};
			my $coor = $orig{$peakpos}{coor};
			my $xpos = $orig{$peakpos}{xpos};
			my $ypos = $orig{$peakpos}{ypos};
			next if defined $used{orig12}{$orig2};
			next if defined $used{orig12b}{$orig1} and $used{orig12b}{$orig1} > $fraction;
			$best{orig12}{$orig2}{peak} = $orig1;
			$best{orig12}{$orig2}{diff} = $diff;
			$best{orig12}{$orig2}{xpos} = $xpos;
			$best{orig12}{$orig2}{ypos} = $ypos;
			$best{orig12}{$orig2}{coor} = $coor;
			$used{orig12}{$orig2} = 1;
			$used{orig12b}{$orig1} ++;
		}
		$fraction = int(@orig2 / @shuf1) + 1;
		foreach my $peakpos (sort {$shuf1{$a}{diff} <=> $shuf1{$b}{diff} || $shuf1{$a}{peak1} <=> $shuf1{$b}{peak1} || $shuf1{$a}{peak2} <=> $shuf1{$b}{peak2}} keys %shuf1) {
			my ($orig2, $shuf1) = split(",", $peakpos);
			my $diff = $shuf1{$peakpos}{diff};
			my $coor = $shuf1{$peakpos}{coor};
			my $xpos = $shuf1{$peakpos}{xpos};
			my $ypos = $shuf1{$peakpos}{ypos};
			next if defined $used{shuf11}{$orig2};
			next if defined $used{shuf11b}{$shuf1} and $used{shuf11b}{$shuf1} > $fraction;
			$best{shuf11}{$orig2}{peak} = $shuf1;
			$best{shuf11}{$orig2}{diff} = $diff;
			$best{shuf11}{$orig2}{xpos} = $xpos;
			$best{shuf11}{$orig2}{ypos} = $ypos;
			$best{shuf11}{$orig2}{coor} = $coor;
			$used{shuf11}{$orig2} = 1;
			$used{shuf11b}{$shuf1} ++;
		}

		$fraction = int(@shuf1 / @orig2) + 1;
		foreach my $peakpos (sort {$shuf1{$a}{diff} <=> $shuf1{$b}{diff} || $shuf1{$a}{peak1} <=> $shuf1{$b}{peak1} || $shuf1{$a}{peak2} <=> $shuf1{$b}{peak2}} keys %shuf1) {
			my ($orig2, $shuf1) = split(",", $peakpos);
			my $diff = $shuf1{$peakpos}{diff};
			my $coor = $shuf1{$peakpos}{coor};
			my $xpos = $shuf1{$peakpos}{xpos};
			my $ypos = $shuf1{$peakpos}{ypos};
			next if defined $used{shuf12}{$shuf1};# or defined $used{shuf12}{$shuf1};
			$best{shuf12}{$shuf1}{peak} = $orig2;
			$best{shuf12}{$shuf1}{diff} = $diff;
			$best{shuf12}{$shuf1}{xpos} = $xpos;
			$best{shuf12}{$shuf1}{ypos} = $ypos;
			$best{shuf12}{$shuf1}{coor} = $coor;
			$used{shuf12}{$shuf1} = 1;
#			print "2: $orig2 -> $orig1\n";
		}
		foreach my $peakpos (sort {$shuf2{$a}{diff} <=> $shuf2{$b}{diff} || $shuf2{$a}{peak1} <=> $shuf2{$b}{peak1} || $shuf2{$a}{peak2} <=> $shuf2{$b}{peak2}} keys %shuf2) {
			my ($orig1, $shuf2) = split(",", $peakpos);
			my $diff = $shuf2{$peakpos}{diff};
			my $coor = $shuf2{$peakpos}{coor};
			my $xpos = $shuf2{$peakpos}{xpos};
			my $ypos = $shuf2{$peakpos}{ypos};
			next if defined $used{shuf21}{$orig1};# or defined $used{2}{$orig2};
			$best{shuf21}{$orig1}{peak} = $shuf2;
			$best{shuf21}{$orig1}{diff} = $diff;
			$best{shuf21}{$orig1}{xpos} = $xpos;
			$best{shuf21}{$orig1}{ypos} = $ypos;
			$best{shuf21}{$orig1}{coor} = $coor;
			$used{shuf21}{$orig1} = 1;
			#$used{2}{$orig2} = 1;
#			print "1: $orig1 -> $orig2\n";
		}
		foreach my $peakpos (sort {$shuf2{$a}{diff} <=> $shuf2{$b}{diff} || $shuf2{$a}{peak1} <=> $shuf2{$b}{peak1} || $shuf2{$a}{peak2} <=> $shuf2{$b}{peak2}} keys %shuf2) {
			my ($orig1, $shuf2) = split(",", $peakpos);
			my $diff = $shuf2{$peakpos}{diff};
			my $coor = $shuf2{$peakpos}{coor};
			my $xpos = $shuf2{$peakpos}{xpos};
			my $ypos = $shuf2{$peakpos}{ypos};
			next if defined $used{shuf22}{$shuf2};# or defined $used{shuf22}{$shuf2};
			$best{shuf22}{$shuf2}{peak} = $orig1;
			$best{shuf22}{$shuf2}{diff} = $diff;
			$best{shuf22}{$shuf2}{xpos} = $xpos;
			$best{shuf22}{$shuf2}{ypos} = $ypos;
			$best{shuf22}{$shuf2}{coor} = $coor;
			$used{shuf22}{$shuf2} = 1;
#			print "2: $orig2 -> $orig1\n";
		}
		my $used1tot = (keys %{$used{orig11}}); my $orig1tot = @orig1;

