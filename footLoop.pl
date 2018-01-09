#!/usr/bin/perl
# Version 160831_Fixed_PrintOutput at the same file (step 8)
use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars   qw($opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f);
my @opts = qw($opt_r $opt_g $opt_i $opt_n $opt_L $opt_x $opt_y $opt_p $opt_q $opt_Z $opt_h $opt_H $opt_F $opt_f);
getopts("r:g:i:n:L:x:y:q:HhZFfp");
BEGIN {
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";

###################
# 0. Check Sanity #
##################

my %opts = ("r" => $opt_r, "g" => $opt_g, "i" => $opt_i, "n" => $opt_n, "L" => $opt_L, "x" => $opt_x, "y" => $opt_y, "p" => $opt_p, "q" => $opt_q, "Z" => $opt_Z, "F" => "NONE", "f" => "NONE", "Z" => "NONE", "p" => "NONE", "L" => $opt_L, "q" => $opt_q);
my %opts2 = ("p" => $opt_p, "Z" => $opt_Z, "F" => $opt_F, "f" => $opt_f);
sanityCheck(\%opts, \%opts2);

###################
# 1. Define Input #
###################
my ($readFile, $genomeFile, $geneIndexesFile, $outDir) = getFullpathAll($opt_r, $opt_g, $opt_i, $opt_n);
my ($readFilename) = getFilename($opt_r, "full");
my $minMapQ  = (not defined($opt_q)) ? 0 : $opt_q;
my $minReadL = (not defined($opt_L)) ? 50 : $opt_L =~ /p$/i ? $opt_L =~ /^(.+)p$/i : $opt_L;
my $bufferL  = (not defined($opt_x)) ? 0 : $opt_x;
my $bufferR  = (not defined($opt_y)) ? 0 : $opt_y;
my $x = defined($opt_x) ? $opt_x : 0;
my $y = defined($opt_y) ? $opt_y : 0;
my $readName = getFilename($readFile, "full");
my $bismarkOpt = (defined $opt_Z) ? "--bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8" : "--bowtie2 --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3";
my $mysam = ($readFile =~ /.f(ast)?q(.gz)?$/) ? $outDir .  "/$readName\_bismark_bt2.sam" : $readFile;
my $uuid = getuuid();
my $date = getDate();

# Make directory
makedir($outDir);

# Make log file
my $logFile = "$outDir/logFile.txt";
open(my $outLog, '>', $logFile);

# Record all options
record_options(\%opts, \%opts2, $outLog);

###################
# 2. Runs Bismark #
###################


# Get gene names from $geneIndexesFile
$geneIndexesFile = get_geneIndexes_fasta($geneIndexesFile, $outDir, $logFile, $outLog);
my ($geneIndexesHash, $geneIndexesFa, $bismark_folder) = parse_geneIndexesFile($geneIndexesFile, $outDir, $outLog);
my %geneIndexes = %{$geneIndexesHash};
make_bismark_index($geneIndexesFa, $bismark_folder, $bismarkOpt, $outLog);

# Run Bismark
run_bismark($readFile, $outDir, $mysam, $opt_F, $outLog);

# Takes sequence of gene and splits into array of individual bases
my $seqFile = $geneIndexesFa;
print $outLog "
!seq=$seqFile
!sam=$mysam
";
my %seq;
print STDERR "\n$YW-------------->$N\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n$YW<--------------$N\n";
print $outLog "\n$YW-------------->$N\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n$YW<--------------$N\n";
open(my $SEQIN, "<", $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
my $fasta = new FAlite($SEQIN);
my %lotsOfC;
while (my $entry = $fasta->nextEntry()) {
	my $gene = uc($entry->def);
	my $seqz = uc($entry->seq);
	$gene =~ s/^>//;
	@{$seq{$gene}{seq}} = split("", $seqz);
	$seq{$gene}{geneL} = @{$seq{$gene}{seq}};
	$seq{$gene}{minReadL} = (defined $opt_L and $opt_L =~ /p$/i) ? int(0.5+$seq{$gene}{geneL} * $minReadL / 100) : $minReadL;
	$seq{$gene}{total} = 0;
	$seq{$gene}{badlength} = 0;
	$seq{$gene}{lowq} = 0;
	$seq{$gene}{used} = 0;
	$seq{$gene}{pos} = 0;
	$seq{$gene}{neg} = 0;
	$seq{$gene}{orig} = $gene;
	print STDERR "\t\tgenez=$gene ($gene) Length=$seq{$gene}{geneL}\n";
	print $outLog "\t\tgenez=$gene ($gene) Length=$seq{$gene}{geneL}\n";
	$lotsOfC{$gene} = get_lotsOfC($seqz);
}
#push(@{$box{$chr}}, "$chr\t$beg\t$end\t$gene\t$val\t$strand$others");
close $SEQIN;
print STDERR "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";
print $outLog "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";

#=OPT_E

my ($samMD5) = getMD5($mysam);
my $origDir = "$outDir/.0_orig_$samMD5/";
check_if_result_exist(["$origDir/.GOOD"], $outLog);
my $checkSam = 1;
my ($mysamName) = getFilename($mysam, "full");
$checkSam = 0 if not -e "$origDir/$mysamName.fixed" and not -e "$origDir/$mysamName.fixed.gz";
makedir($origDir);
if (-e "$origDir/$mysamName.fixed.gz") {
	my ($samLineCount2) = linecount("$origDir/$mysamName.fixed.gz");
	my ($samLineCount1) = linecount($mysam);
	my $nine = $samLineCount1 * 0.9;
	$checkSam = 0 if $nine > $samLineCount2;
}
elsif (-e "$origDir/$mysamName.fixed") {
	my ($samLineCount2) = linecount("$origDir/$mysamName.fixed");
	my ($samLineCount1) = linecount($mysam);
	$checkSam = 0 if $samLineCount1 - 500 > $samLineCount2;
	print "gzip $origDir/$mysamName.fixed\n";
	print "md5sum $origDir/$mysamName.fixed.gz > $origDir/.$mysamName.fixed.gz.md5\n";
	system("gzip $origDir/$mysamName.fixed") == 0 or die "Failed to gzip $origDir/$mysamName.fixed: $!\n";
	system("md5sum $origDir/$mysamName.fixed.gz > $origDir/.$mysamName.fixed.gz.md5") == 0 or die "Failed to md5sum $origDir/$mysamName.fixed > $origDir/.$mysamName.fixed.gz.md5: $!\n";
}
if ($checkSam == 0) {
	print STDERR "line 122 of footLoop.pl shouldnt be accesse\n";
	print $outLog "line 122 of footLoop.pl shouldnt be accesse\n" and exit 1;
	print STDERR "footLoop_2_sam_to_peak.pl -f $outDir -o $origDir\n";
	print $outLog "footLoop_2_sam_to_peak.pl -f $outDir -o $origDir\n";
	system("footLoop_2_sam_to_peak.pl -f $outDir -o $origDir") == 0 or die "Failed to run footLoop_2_sam_to_peak.pl -f $outDir: $!\n";
	system("gzip $origDir/$mysamName.fixed.gz") == 0 or die "Failed to gzip $origDir/$mysamName.fixed.gz: $!\n";
}

# File description, open files, etc
my $samFile = "$mysam";
open(my $notused, ">", "$outDir/.$readFilename.notused") or die "Cannot open $outDir/.$readFilename.notused: $!\n";
print STDERR "\n$YW-------------->$N\n${YW}3. Parsing sam file $CY$samFile$N and getting only high quality reads\n$YW<--------------$N\n";
print $outLog "\n$YW-------------->$N\n${YW}3. Parsing sam file $CY$samFile$N and getting only high quality reads\n$YW<--------------$N\n";
open(my $sam, $samFile) or print $outLog "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!" and exit 1;

## Some stats
my $linecount = 0; my %count; ($count{total}, $count{used}, $count{diffgene}, $count{lowq}, $count{badlength}) = (0,0,0,0,0);
my %readz;

#####################
# A. Parse Sam File #
#####################
# Loops through each read and writes high quality reads into two separate files <gene>Positive.txt and <gene>Negative.txt

while(my $line = <$sam>) {
	$linecount ++;
	chomp($line);
	
	#####################
	# 1. Parse sam line #
	#####################
	my @arr = split("\t", $line);
	print STDERR "\t$YW$samFile$N: Done $linecount\n" if $linecount % 5000 == 0;
	print $outLog "\t$YW$samFile$N: Done $linecount\n" if $linecount % 5000 == 0;

	# a. Total column must be 14 or skipped (e.g. sam header)
	print $outLog "$LGN$linecount$N: SAM header detected (#column < 14). LINE:\n\t$line\n" and next if @arr < 14;
	
	my ($eval, $evalPrint) = myeval(\@arr);
	my ($read, $readStrand, $chrom, $readPos, $mapQ, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $tags, $junk4, $readMethCall) = @arr;
	my $cigarLen = parse_cigar($cigar, "len");
	my $seqsLen  = length($seqs);
	my $gene     = uc($chrom);

	print $outLog "$LGN$linecount$N: read $LCY$read$N is mapped to gene $LGN$gene$N is not in geneIndexes! LINE:\n\t$line\n" if not defined $seq{$gene};
	next if not defined($seq{$gene});

	my $cigarPerc = int($cigarLen / $seqsLen * 10000)/100;
	next if $cigarPerc < 75;
	
	# (EXTRA) This below is to show the user an example of parsed line and tells user if it's parsed every 20k read.
	my ($readname) = length($read) >= 20 ? $read =~ /(.{20})$/ : $read; $readname = "..." . $readname  if length($read) >= 20;
	$count{total} ++;
	$seq{$gene}{total} ++;
	print $outLog "\tExample at read $count{total}: name=$CY$readname$N\tstrand=$CY$readStrand$N\tchr/gene=$CY$chrom$N\tpos=$CY$readPos$N\tmapQ=$CY$mapQ$N\n\n" if $count{total} == 1;
	print $outLog "\tDone $GN$count{total}$N\n" if $count{total} % 20000 == 0;
	
	# b. Filter out reads with mapping quality=$mapQ less than threshold quality=$opt_q
	if($mapQ < $opt_q)
	{
		print $notused "\t$CY$readname$N quality ($CY$mapQ$N) is less than $CY$opt_q$N\n";
		$count{lowq} ++;
		$seq{$gene}{lowq} ++;
		next;
	}
	
	# c. Filter out reads with read length=length($seqs) less than  the proper length, or at the proper location in genome(accounting for indexing)
	elsif(length($seqs) < $seq{$gene}{minReadL}) # buffer is obsolete
	{
		# (EXTRA) This below is just for statistics
		$count{badlength} ++;
		$seq{$gene}{badlength} ++;
		die "READNAME undef\n" if not defined($readname);
		die "fields9 undef\n" if not defined($seqs);
		die "length of gene $chrom undef\n" if not defined($seq{$gene}{geneL});
		
		print $notused "\t$CY$readname$N length of seq ($CY" . length($seqs). "$N) is less than 500bp (length of original sequence is ($CY" . $seq{$gene}{geneL} . "$N)!\n";
		next;
	}
	
	#counts number of CT conversions (takes into account whether or not the user wants to include Cs in CpG context)
	#my $CT = ($readMethCall =~ tr/xhu/xhu/);
	#if($opt_c)
	#{
	#	$CT = ($readMethCall =~ tr/zxhu/zxhu/);
	#}
	
	#writes positive reads into <gene>Positive.txt and negative reads into <gene>Negative.txt
	if($readStrand == 0 || $readStrand == 16)
	{
		#my $to_be_printed = "$CT\t$line\n";
		#$readStrand == 0 ? print $positiveReads $to_be_printed : print $negativeReads $to_be_printed;
		$seq{$gene}{read}{$read} = $readStrand;
	}
}


# print from insamfix
# Open the .orig files
makedir("$origDir") if not -d "$origDir";
foreach my $genez (sort keys %seq) {
	my $outTXTFilePos  = "$origDir/$genez\_Pos.orig";
	my $outTXTFileNeg  = "$origDir/$genez\_Neg.orig";
	my $outTXTFileUnk  = "$origDir/$genez\_Unk.orig";
#	my $outTXTFilePosR = $outDir . "0_orig/$genez\_RevPos.orig";
#	my $outTXTFileNegR = $outDir . "0_orig/$genez\_RevNeg.orig";
	open ($seq{$genez}{outTXTPos}, ">", $outTXTFilePos);
	open ($seq{$genez}{outTXTNeg}, ">", $outTXTFileNeg);
	open ($seq{$genez}{outTXTUnk}, ">", $outTXTFileUnk);
#	open ($seq{$genez}{outTXTRevPos}, ">", $outTXTRevFilePos);
#	open ($seq{$genez}{outTXTRevNeg}, ">", $outTXTRevFileNeg);
}

my $skipped = 0; my ($passedFilterP, $passedFilterN) = (0,0);
open (my $inSamFix, "zcat $origDir/$mysamName.fixed.gz|") or print $outLog "Failed to open $origDir/$mysamName.fixed.gz: $!\n" and exit 1;
while (my $line = <$inSamFix>) {
	chomp($line);
	my ($read,$type, $oldStrand, $strand, $genez, $pos, $info) = split("\t", $line);
	my ($CT0, $CC0, $GA0, $GG0, $CT1, $CC1, $GA1, $GG1) = split(",", $info);
	if (not defined $seq{$genez} or not defined $info) {
		print $outLog "\tERROR in $LCY$origDir/$mysamName.fixed$N: gene=$genez but \$seq{\$genez} is not defined!line=\n$line\n\n" if not defined $seq{$genez};
		print $outLog "\tERROR in $LCY$origDir/$mysamName.fixed$N: gene=$genez and seq genez is $seq{$genez} but info is not defined!line=\n$line\n\n" if defined $seq{$genez} and not defined $info;
		$skipped ++;
		next;
	}
	if (not defined $seq{$genez}{read}{$read}) {
		next;
	}
#	my $oldStrand = $seq{$genez}{read}{$read};
	if ($type eq "6_BOTH") {
		print {$seq{$genez}{outTXTUnk}} "$read\tBP\t$pos\n" if $strand eq 0;
		print {$seq{$genez}{outTXTUnk}} "$read\tBN\t$pos\n" if $strand eq 16;
		$seq{$genez}{unkpos} ++ if $strand == 0;
		$seq{$genez}{unkneg} ++ if $strand == 16;
		$passedFilterP ++ if $strand == 0;
		$passedFilterN ++ if $strand == 16;
	}
	else {
		print {$seq{$genez}{outTXTNeg}} "$read\tFN\t$pos\n" if $strand eq 16;
		print {$seq{$genez}{outTXTPos}} "$read\tFP\t$pos\n" if $strand eq 0;
		$seq{$genez}{pos} ++ if $strand == 0;
		$seq{$genez}{neg} ++ if $strand == 16;
		$passedFilterP ++ if $strand == 0;
		$passedFilterN ++ if $strand == 16;
	}
	$count{used} ++;
	$seq{$genez}{used} ++;
	$seq{$genez}{posneg} ++ if $strand eq 16 and $oldStrand eq 0;
	$seq{$genez}{negpos} ++ if $strand eq 0 and $oldStrand eq 16;
}
close $inSamFix;
print $outLog "\t$LCY$origDir/$mysamName.fixed$N: skipped = $LGN$skipped$N\n";

###########
#=comment
#      my $CTval = $CT eq " "       ? "-"  : $CT eq "."       ? $ref     : $CT =~ /^[Nn]$/ ? $CT :
#                  $CT =~ /^[Z]$/   ? "Z"  : $CT =~ /^[z]$/   ? "z"      : $CT eq "-" ? "-" :
#                  $CT =~ /^[XHU]$/ ? $ref : $CT =~ /^[xhu]$/ ? lc($ref) : die "Cannot determine CT=$CT\n";
#
#      my $GAval = $GA eq " "       ? "-"  : $GA eq "."       ? $ref     : $GA =~ /^[Nn]$/ ? $GA :
#                  $GA =~ /^[Z]$/   ? "Z"  : $GA =~ /^[z]$/   ? "z"      : $GA eq "-" ? "-" :
#                  $GA =~ /^[XHU]$/ ? $ref : $GA =~ /^[xhu]$/ ? lc($ref) : die "Cannot determine GA=$GA\n";
#=cut

print $outLog "
Reads that passed filters:
Positive: $passedFilterP
Negative: $passedFilterN
Total   : $count{total};

Per Gene:
";

my $zero = "";
foreach my $gene (sort keys %seq) {
	my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
	foreach my $key (@key) {
		$seq{$gene}{$key} = 0 if not defined $seq{$gene}{$key};
	}
}
foreach my $gene (sort {$seq{$b}{total} <=> $seq{$a}{total}} keys %seq) {
	my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
	my $outTXTFilePos  = "$origDir/$gene\_Pos.orig"; system("/bin/rm $outTXTFilePos") if -e $outTXTFilePos and -s $outTXTFilePos == 0;
	my $outTXTFileNeg  = "$origDir/$gene\_Neg.orig"; system("/bin/rm $outTXTFileNeg") if -e $outTXTFileNeg and -s $outTXTFileNeg == 0;
	my $outTXTFileUnk  = "$origDir/$gene\_Unk.orig"; system("/bin/rm $outTXTFileUnk") if -e $outTXTFileUnk and -s $outTXTFileUnk == 0;

	$zero .= "$gene ($seq{$gene}{total})\n" and next if $seq{$gene}{total} <= 10;
	my $gene2 = $seq{$gene}{orig};
	print $outLog "
- $gene (original name = $gene2):
Positive    = $seq{$gene}{pos} (neg->pos = $seq{$gene}{negpos})
Negative    = $seq{$gene}{neg} (pos->neg = $seq{$gene}{posneg})
UnkPos      = $seq{$gene}{unkpos}
UnkNeg      = $seq{$gene}{unkneg}
Used        = $seq{$gene}{used}
Total       = $seq{$gene}{total}
Too Short   = $seq{$gene}{badlength}
Low Quality = $seq{$gene}{lowq}
";

	print STDERR "
- $gene (original name = $gene2):
Positive    = $seq{$gene}{pos} (neg->pos = $seq{$gene}{negpos})
Negative    = $seq{$gene}{neg} (pos->neg = $seq{$gene}{posneg})
UnkPos      = $seq{$gene}{unkpos}
UnkNeg      = $seq{$gene}{unkneg}
Used        = $seq{$gene}{used}
Total       = $seq{$gene}{total}
Too Short   = $seq{$gene}{badlength}
Low Quality = $seq{$gene}{lowq}
";
}
$zero = $zero eq "" ? "(None)\n" : "\n$zero\n";
print STDERR "- Genes that have <= 10 reads:$zero\n";
print STDERR "\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";
print $outLog "\n\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";

print STDERR "\n\nOutput: ${LGN}$origDir$N\n";
print $outLog "\n\nOutput: ${LGN}$origDir$N\n";
##########################################

foreach my $gene (sort keys %seq) {
	my @seq = @{$seq{$gene}{seq}};
	my $finalPos = $origDir . "/$gene\_Pos.orig";
	my $finalNeg = $origDir . "/$gene\_Neg.orig";

#=PART6
	my $lotsOfC = defined $lotsOfC{$gene} ? $lotsOfC{$gene} : "";
	my $ucgene = uc($gene);
	print "footLoop_addition.pl $finalPos $ucgene $lotsOfC\n";
	print "footLoop_addition.pl $finalNeg $ucgene $lotsOfC\n";
	#system("footLoop_addition.pl $finalPos $ucgene \"$lotsOfC\"") == 0 or print "Cannot do footLoop_addition.pl $ucgene $finalPos: $!\n";
}
die;
=comment
	print $outLog "\tfootLoop_addition.pl $finalPos $ucgene \"$lotsOfC\"\n";
	my ($finalPosLine) = `wc -l $finalPos` =~ /^(\d+) /;
	my ($finalNegLine) = `wc -l $finalNeg` =~ /^(\d+) /;
	
	print $outLog "\n\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";
	print STDERR "\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";


	
	#makes heatmaps for positive and negative strand
	
	#0= not converted (lightyellow) # same
	#1= converted (green)   # same
	#2= non-C (white)       # same
	#3= non-converted CpG (black) # is non converted
	#4= converted CpG  (blue) # is converted
	#5= converted CpG in R-loop (purple) # conv CPG
	#6= no data (grey)
	#9= converted C in R-loop (red) # not exist
	my $threshold = $opt_t * 100;
	my $finalPosPDF = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . ".pdf";
	my $finalNegPDF = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . ".pdf";
	$finalPosPDF = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "CG.pdf" if($opt_c);
	$finalNegPDF = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "CG.pdf" if($opt_c);
	my $finalPosPDF_NOPEAK = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "_NOPEAK.pdf";
	my $finalNegPDF_NOPEAK = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "_NOPEAK.pdf";
	$finalPosPDF_NOPEAK = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "_NOPEAK_CG.pdf" if($opt_c);
	$finalNegPDF_NOPEAK = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "_NOPEAK_CG.pdf" if($opt_c);
	my $finalPosPNG = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . ".png";
	my $finalNegPNG = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . ".png";
	$finalPosPNG = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "CG.png" if($opt_c);
	$finalNegPNG = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "CG.png" if($opt_c);
	my $finalPosPNG_NOPEAK = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "_NOPEAK.png";
	my $finalNegPNG_NOPEAK = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "_NOPEAK.png";
	$finalPosPNG_NOPEAK = $outDir . "$outDirname\_$gene\_Pos" . ($opt_t*100) . "_NOPEAK_CG.png" if($opt_c);
	$finalNegPNG_NOPEAK = $outDir . "$outDirname\_$gene\_Neg" . ($opt_t*100) . "_NOPEAK_CG.png" if($opt_c);

	my $Rscript = "$outDir/$outDirname\_$gene\_MakeHeatmap.R";
	open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";

	# Breaks and color for R script
	my $breaks = "c(-0.5"	; my $colors = "c("				 ;
	$breaks .= ",0.5"			; $colors .= "\"cornsilk\""	 ; # 0 = not converted
	$breaks .= ",1.5"			; $colors .= ",\"green4\""		 ; # 1 = converted C
	$breaks .= ",2.5"			; $colors .= ",\"white\""		 ; # 2 = A T or G (non C)
	if (defined $opt_c) {
		$breaks .= ",3.5"		; $colors .= ",\"peachpuff\""  ; # 3 = Non converted CpG
		$breaks .= ",4.5"		; $colors .= ",\"seagreen4\""  ; # 4 = Converted CpG
		$breaks .= ",5.5"		; $colors .= ",\"maroon4\""    ; # 5 = PEAK Converted CpG
	}
	$breaks .= ",6.5"			; $colors .= ",\"grey\""		 ; # 6 = No data
	$breaks .= ",9.5"			; $colors .= ",\"red4\""		 ; # 9 = PEAK Converted C
	# For nucleotide
	#$breaks .= ",10.5"		; $colors .= ",\"seagreen4\""  ; # 10 = Nucleotide A
	#$breaks .= ",11.5"		; $colors .= ",\"royalblue3\"" ; # 11 = Nucleotide C
	#$breaks .= ",12.5"		; $colors .= ",\"saddlebrown\""; # 12 = Nucleotide T
	#$breaks .= ",13.5"		; $colors .= ",\"red2\""       ; # 13 = Nucleotide G
	$breaks .= ")"    		; $colors .= ")"               ;

	#Function to turn @seq into number and array for R
	my $seqR = "seq = c(";
	for (my $s = 0; $s < @seq; $s++) {
		my $nuc = $seq[$s];
		#$nuc = $nuc =~ /A/i ? 11 : $nuc =~ /C/i ? 12 : $nuc =~ /G/i ? 13 : $nuc =~ /T/i ? 14 : 15;
		$nuc = $nuc =~ /A/i ? "seagreen4" : $nuc =~ /C/i ? "blue4" : $nuc =~ /G/i ? "saddlebrown" : $nuc =~ /T/i ? "red4" : "grey";

		$seqR .= "\"$nuc\"";
		$seqR .= $s == @seq - 1 ? ")" : ",";
	}
	
	my $labelz = "";
	if (defined $opt_T) {
		$labelz .= "p = p + geom_text(aes(x=min(dm\$x),label=id),hjust=0)";
	}
	my $genez = uc($gene);
	my $RBox = ""; my $RBoxPlot = "";
	if (defined $box) {
		if (defined $box->{$genez}) {
			my %RBox;
			foreach my $lines (sort @{$box->{$genez}}) {
				my ($chr, $beg, $end, $geneBox) = split("\t", $lines);
				#my $beg2 = $geneIndexes{$genez};
				#if ($beg > $beg2) {
				#	$beg -= $beg2;
				#	$end -= $beg2;
				#}
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
		files = c(\"$finalPos\",\"$finalNeg\",\"$finalPos_NOPEAK\",\"$finalNeg_NOPEAK\")
		PDFout = c(\"$finalPosPDF\",\"$finalNegPDF\",\"$finalPosPDF_NOPEAK\",\"$finalNegPDF_NOPEAK\")
		PNGout = c(\"$finalPosPNG\",\"$finalNegPNG\",\"$finalPosPNG_NOPEAK\",\"$finalNegPNG_NOPEAK\")
		PEAKS = paste(files,\".PEAKS\",sep=\"\")

		info = file.info(files)
		files = rownames(info[info\$size != 0, ])
		for (i in 1:length(files)) 
		{
			if (file.exists(files[i])) 
			{
				myclust = TRUE
				if(length(grep(\"_NOPEAK\", files[i])) != 0) {
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
}
print STDERR "\n${LPR}If there is not PDF made from step (8) then it's due to too low number of read/peak$N\n\n";
=cut

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

sub run_bismark {
	my ($readFile, $outDir, $mysam, $force, $outLog) = @_;
	if (defined $force or not -e $mysam) {
		my $text = "\n$YW-------------->$N\n${YW}2. Running bismark$N $bismarkOpt $bismark_folder $readFile\n";
		print STDERR $text; print $outLog $text;
		if ($opt_p) {
			my $outFolder = $outDir . "/.bismark_paralel/";
			if (-d $outFolder) {
				print "Removing all in $outFolder\n";
				system("/bin/rm $outFolder/*");
			}
			makedir($outFolder);
			print $outLog "\t  Splitting $CY$readFile$N by 1000 sequences!\n\n####### SPLITRESULT LOG ########";	
			my $splitresult = `SplitFastq.pl -i $readFile -o $outFolder -n 1000`;
			print $outLog "$splitresult\n";	
			print "####### SPLITRESULT LOG #######\n\n\t  Running bismark in paralel!\n";
			my $result = system("run_script_in_paralel2.pl -v \"srun -p high --mem 8000 bismark -o $outDir/.bismark_paralel/ $bismarkOpt $bismark_folder FILENAME >> FILENAME.bismark.log 2>&1\" $outFolder .part 20");
			my @partSam = <$outFolder/*.part_bismark*.sam>; my $totalPartSam = @partSam;
			print $outLog "\t  All part.sam has been made (total = $totalPartSam). Now making $CY$mysam$N and then removing the part sam\n";
			for (my $p = 0; $p < @partSam; $p++) {
				my $partSam = $partSam[$p];
				print "\t\tPutting $partSam into $mysam and removing it!\n";
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p == 0;
				system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p != 0;
				print "\t- Removing $CY$partSam$N: /bin/rm $partSam\n";
				system("/bin/rm $partSam") == 0 or die "Failed to /bin/rm $partSam: $!\n";
			}
		}
		else {
			my $result = system("bismark -o $outDir $bismarkOpt $bismark_folder $readFile > $outDir/.bismark_log 2>&1");
			if ($result != 0) {
				print STDERR  "\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_folder\n";
				print $outLog "\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_folder\n";
				system("bismark_genome_preparation --bowtie2 $bismark_folder") == 0 or die "Failed to run bismark genome preparation: $!\n";
				system("bismark $bismarkOpt $bismark_folder $readFile") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
			}
			print $outLog "\t${GN}SUCCESS$N: Output $mysam\n";
			print STDERR "\t${GN}SUCCESS$N: Output $mysam\n";
		}
	}
	else {
		print $outLog "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
		print STDERR "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
	}
}

sub make_bismark_index {
	my ($geneIndexFa, $bismark_folder, $bismarkOpt, $outLog) = @_;
	print STDERR  "\n$YW-------------->$N\n${YW}1b. Running bismark_genome_preparation$N --bowtie2 $bismark_folder$N\n$YW<--------------$N\n";
	print $outLog "\n$YW-------------->$N\n${YW}1b. Running bismark_genome_preparation$N --bowtie2 $bismark_folder$N\n$YW<--------------$N\n";
	my ($check, $md5sum, $md5sum2) = (0);
	my $bismark_folder_exist = (-d "$bismark_folder/Bisulfite_Genome/" and -e "$bismark_folder/Bisulfite_Genome/.md5sum.md5") ? 1 : 0;
	if ($bismark_folder_exist == 1) {
		print STDERR "\tOlder bismark folder $LCY$bismark_folder/BIsulfite_Genome/$N exist! Checking MD5 if they're the same as current fasta file.\n";
		print $outLog "\tOlder bismark folder $LCY$bismark_folder/BIsulfite_Genome/$N exist! Checking MD5 if they're the same as current fasta file.\n";
		($md5sum)  = getMD5("$bismark_folder/Bisulfite_Genome/.md5sum.md5");
		($md5sum2) = getMD5($geneIndexesFa);
		$bismark_folder_exist = 0 if $md5sum ne $md5sum2;
	}
	if ($bismark_folder_exist == 0) {
		print STDERR "\tOlder bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if defined $md5sum;
		print $outLog "\tOlder bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if defined $md5sum;
		print $outLog "\tbismark_genome_preparation --bowtie2 $bismark_folder && md5sum $geneIndexesFa > $bismark_folder/Bisulfite_Genome/.md5sum.md5\n";
		system("bismark_genome_preparation --bowtie2 $bismark_folder > $bismark_folder/LOG.txt 2>&1 && md5sum $geneIndexesFa > $bismark_folder/Bisulfite_Genome/.md5sum.md5") == 0 or die "Failed to run bismark genome preparation: $!\n";
	}
	else { 
		print STDERR "\tOlder bisulfite genome have the ${LGN}same$N MD5sum! (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if $md5sum eq $md5sum2;
		print $outLog "\tOlder bisulfite genome have the ${LGN}same$N MD5sum! (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if $md5sum eq $md5sum2;
		print STDERR "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist and is used!\n";
		print $outLog "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist and is used!\n";
	}
}

sub parse_geneIndexesFile {
	my ($geneIndexesFile, $outDir, $outLog) = @_;
	my $geneIndexes;
	open (my $geneIndexIn, "<", $geneIndexesFile) or die "Cannot read from $geneIndexesFile: $!\n";
	while (my $line = <$geneIndexIn>) {
		chomp($line);
		my ($chr, $beg, $end, $gene) = split("\t", $line);
		$gene = uc($gene);
		print $outLog "\tgeneIndexesFile=$geneIndexesFile,gene=$gene,beg=$beg,end=$end\n";
		$geneIndexes->{$gene} = $beg;
	}
	close $geneIndexIn;
	my $geneIndexesName = getFilename($geneIndexesFile);
	my $geneIndexesFaTemp = $outDir . "/.geneIndexes/$geneIndexesName.fa";
	makedir($geneIndexesFaTemp, 1);

	my $text = "\n$YW-------------->$N\n${YW}1a. Getting fasta sequence from $geneIndexesFile into $CY$geneIndexesFaTemp$N\n$YW<--------------$N\n";
	print STDERR $text; print $outLog $text; 
	print $outLog "\t- Running ${YW}bedtools getfasta$N -fi $opt_g -bed $geneIndexesFile -fo $geneIndexesFaTemp -name\n";

	system("fastaFromBed -fi $opt_g -bed $geneIndexesFile -fo $geneIndexesFaTemp -name") == 0 ? print $outLog "\t${GN}SUCCESS$N: Output: $CY$geneIndexesFaTemp$N\n" : die "Failed to run bedtools: $!\n";
	$geneIndexesFaTemp = uppercaseFasta($geneIndexesFaTemp, $outLog);
	my ($geneIndexesFaMD5, $temp, $geneIndexesFaTempMD5File)  = getMD5($geneIndexesFaTemp);
	die "GENE=$geneIndexesFaTemp, MD5=$geneIndexesFaMD5\n" if $geneIndexesFaMD5 eq "NA";
	my $geneIndexesFa = $outDir . "/.geneIndexes/$geneIndexesFaMD5/$geneIndexesName.fa";
	my $geneIndexesFaMD5File = $outDir . "/.geneIndexes/$geneIndexesFaMD5/.$geneIndexesName.fa.md5";
	if (not -d $geneIndexesFaMD5 or (-d $geneIndexesFaMD5 and not -e $geneIndexesFa)) {
		makedir("$outDir/.geneIndexes/$geneIndexesFaMD5");
		system("/bin/mv $geneIndexesFaTemp $outDir/.geneIndexes/$geneIndexesFaMD5/") == 0 or print $outLog "Failed to mv $geneIndexesFaTemp $outDir/.geneIndexes/$geneIndexesFaMD5/: $!\n" and exit 1;
		system("/bin/mv $geneIndexesFaTempMD5File $outDir/.geneIndexes/$geneIndexesFaMD5/") == 0 or print $outLog "Failed to mv $geneIndexesFaTemp $outDir/.geneIndexes/$geneIndexesFaMD5/: $!\n" and exit 1;
	}
	print $outLog "\nERROR: geneIndexesHash is not defined!\n" and exit 1 if not defined $geneIndexes;
	print $outLog "\nERROR: geneIndexesFa is not defined!\n" and exit 1 if not defined $geneIndexesFa;
	return ($geneIndexes, $geneIndexesFa, "$outDir/.geneIndexes/$geneIndexesFaMD5/");
}

sub uppercaseFasta {
	my ($fastaFile, $outLog) = @_;
	open (my $in, "<", $fastaFile) or print $outLog "Failed to open $fastaFile: $!\n" and exit 1;
	open (my $out, ">", "$fastaFile.out") or print $outLog "Failed to write to $fastaFile.out: $!\n" and exit 1;
	my $fasta = new FAlite($in);
	while (my $entry = $fasta->nextEntry()) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$seq =~ tr/a-z/A-Z/;
		print $out "$def\n$seq\n";
	}
	close $in;
	close $out;
	system("/bin/mv $fastaFile.out $fastaFile") == 0 or print $outLog "Failed to /bin/mv $fastaFile.out $fastaFile: $!\n" and exit 1;
	return("$fastaFile");
}
sub get_geneIndexes_fasta {
	my ($geneIndexesFile, $outDir, $logFile, $outLog) = @_;
	my $geneIndexesName = getFilename($geneIndexesFile);
	#print "GENE indexes=$geneIndexesFile, name = $geneIndexesName\n";
	if ($x eq 0 and $y eq 0) {
		system("/bin/cp $geneIndexesFile $outDir/$geneIndexesName\_$x\_$y\_bp.bed") == 0 or print $outLog "Failed to copy $geneIndexesFile to $outDir/$geneIndexesName\_$x\_$y\_bp.bed: $!\n" and exit 1;
	}
	else {
		system("bedtools_bed_change.pl -m -x $x -y $y -i $geneIndexesFile -o $outDir/$geneIndexesName\_$x\_$y\_bp.bed >> $logFile 2>&1") == 0 or print $outLog "Failed to get (beg $x, end $y) bp of $geneIndexesFile!\n" and exit 1;
	}
	$geneIndexesFile = "$outDir/$geneIndexesName\_$x\_$y\_bp.bed";
	return($geneIndexesFile);
}

#=SUB convert_seq {

sub sanityCheck {

my ($opts) = @_;

my $usageshort = "\n$LGN----------------------------$N\n
Usage: $YW$0$N [options..]
\t$CY-r$N read.fq
\t$LPR-n$N output_dir
\t$LGN-g$N genome.fa
\t$YW-i$N geneIndexes.bed
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

my $checkopts = 0;
foreach my $key (sort keys %{$opts}) {
	$checkopts ++ if defined $opts->{$key};
}
die $usageshort if $checkopts == 0;

#Usage: $YW$0$N ${GN}[options]$N [-r $CY<read.fq>$N OR -S $CY<read.bismark.sam>$N -n $PR<output folder>$N -g $GN<genome.fa>$N -i $YW<geneCoordinates.bed>$N
my $usage = "
${LGN}Options [default]:$N
-x $LGN<Left pad [0]>$N -y $LGN<Right pad [0]>$N -t $LGN<threshold [0.65]>$N -l $LGN<min peak length [100]>$N -L $LGN<min READ length [500]>$N -q $LGN<min map quality [0]>$N

Example:
$YW$0$N\n\t$CY-r$N example/pacbio12ccs.fq \n\t$PR-n$N myoutput\n\t$GN-g$N example/hg19.fa.fa\n\t$YW-i$N example/pacbio12index.bed\n\t-x-100\n\t-y 100\n\t-l 100\n\t-L 1000\n\t-t 0.5

$LGN-Z$N: Make mapping parameters not stringent. 
Bismark uses bowtie2 mapping parameters:
- stringent (default) formula: --rdg 5,3 --rfg 5,3 --score_min L,0,-0.3
- not stringent (-Z)  formula: --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8

-H: Longer explanation about bismark parameters above.


(..continued; Use -H to display longer help!)

";

my $usage_long = "$usage

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

die "\n$usage\n" if defined $opt_h;
die "\n$usage_long\n" if defined $opt_H;
my $read0 = defined($opt_r) ? $opt_r : "FALSE";

die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -r <read.fq> not defined\n\n############################################\n\n" if not defined($opt_r);
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -r $opt_r DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if defined $opt_r and not -e $opt_r;
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -n <output directory> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_n);
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -i <geneindex.bed> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_i);
#die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -l <min peak length in bp> must be positive integer!\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if defined($opt_l) and $opt_l !~ /^\d+$/;
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -L <min read length in bp> must be positive integer!\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if defined($opt_L) and $opt_L !~ /p$/ and $opt_L !~ /^\d+$/;
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -g <ref_genome.fa [hg19.fa]> not defined\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not defined($opt_g);
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -i $opt_i DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not -e ($opt_i);
die "\n${LRD}########## ${N}FATAL ERROR${LRD} ##########\n\nREASON:$N -g $opt_g DOES NOT EXIST\n\n############################################\n\nDo $YW$0$N $LGN-h$N to see usage info\n\n" if not -e ($opt_g);

	if (not -d "$footLoopDir/.sortTMP/") {
		system("mkdir $footLoopDir/.sortTMP/") == 0 or die "Failed to make directory $footLoopDir/.sortTMP/: $!\n";
	}

	my ($PositiveFinal, $NegativeFinal, $MethylationPos, $methylationNeg);

	my $outDir = myFootLib::getFullpath("$opt_n") . "/";
	my @dir = split("/", $outDir);
	my $outDirname = $dir[@dir-1];
	   $outDirname = $dir[@dir-2] if not defined $outDirname or (defined $outDirname and $outDirname =~ /^[\s]*$/);

}

sub record_options {
	my ($opts, $opts2, $outLog) = @_;
	my $optPrint = "$0";
	foreach my $opt (sort keys %{$opts}) {
		if (defined $opts->{$opt} and $opts->{$opt} eq "NONE") {
			$optPrint .= " -$opt" if defined $opts2->{$opt};
		}
		elsif (defined $opts->{$opt} and $opts->{$opt} ne "NONE") {
			$optPrint .= " -$opt $opts->{$opt}";
		}
	}
	my  $param = "
${YW}Initializing...$N
	
Date       : $date
Run ID     : $uuid
Run script : $optPrint
	
${YW}Input Parameters$N
1. -r ${CY}Read/SAM$N  :$readFile
2. -g ${CY}Genome$N    :$genomeFile
3. -n ${CY}OutDir$N    :$outDir
4. -i ${CY}Index$N     :$geneIndexesFile
5. -L ${CY}MinRdLen$N  :$minReadL
6. -q ${CY}MinMapQ$N   :$opt_q
7. -Z ${CY}Bismark$N   :$bismarkOpt

";
	print $outLog $param;
	print STDERR $param;
}

sub check_if_result_exist {
	my ($resFiles, $outLog) = @_;
	return if not defined $resFiles;
	my @resFiles = @{$resFiles};	
	return if @resFiles == 0;
	if (defined $opt_n and -d $opt_n and not defined $opt_f) {
		if (@{$resFiles} != 0) {
			my $resFileCheck = 0;
			foreach my $resFile (@{$resFiles}) {
				$resFileCheck = 1 if -e $resFile and -s $resFile > 0;
			}
			print $outLog "-d $LGN$opt_n$N $resFiles->[0] already exist!\n" if $resFileCheck == 1;
			exit 0 if $resFileCheck == 1;
		}
		else {
			print $outLog "End of run file does not exist!\n";
		}
	}
	return;
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
			print $outLog "\t- Doing gene $CY$gene$N ($fileLast)\n";
			print $outLog "$fileLast doesn't exist!\n" and next if not -e $fileLast;
			print $outLog "$fileLast has no line!\n" and next if `wc -l < $fileLast` == 0;
			open (my $fileLastIn, "<", $fileLast) or print $outLog "Cannot open $fileLast\n" and die "Cannot open $fileLast: $!\n";
			my ($fileLastCount) = `wc -l $fileLast` =~ /^(\d+) /;
			
			my %peak;
			my $start = 1;
			my $end = $minPeakLength;
			my $linecount = 0;
			my $LENGTH = 0; my $NAME = "NA";
			while (my $lineFinal = <$fileLastIn>)
			{
				$linecount ++;
				chomp($lineFinal);
				print $outLog "\tDone $GN$linecount$N\n" if $linecount % 500 eq 0;
				print "\t$fileLast: Done $GN$linecount$N / $fileLastCount\n" if $linecount % 500 eq 0;
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
			print $outLog "$finalPos\tlength=$maxAllLength\tpeak=$peakcount\tnopeak=$nopeak\ttotal=$total\n" if $strand == 0;
			print $outLog "$finalNeg\tlength=$maxAllLength\tpeak=$peakcount\tnopeak=$nopeak\ttotal=$total\n" if $strand == 1;
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
	print $outLog "\t${BR}2a. Parsing in exon intron data from $CY$exonFile$N:\n";
	foreach my $gene (sort keys %seq) {
		next if -e "$outFolder/exon/$gene.exon"; # DELETE THIS
		my ($genepos) = `grep -i -P \"\\t$gene\\t\" $geneIndexesFile`; chomp($genepos);
		my $length_seq = $seq{$gene}{geneL};
		print $outLog "\n$LRD!!!$N\tWARNING: No sequence for gene $CY$gene$N is found in -s $CY$seqFile$N!\n\n" if $length_seq == 0;
		myFootLib::parseExon($exonFile, $genepos, $gene, $outFolder, $length_seq);
		print STDERR "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
		print $outLog "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
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
