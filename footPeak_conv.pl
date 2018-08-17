#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_n $opt_i $opt_S $opt_G);
getopts("vn:i:SG:");

#########
# BEGIN #
#########

BEGIN {
	my ($bedtools) = `bedtools --version`;
	my ($bowtie2) = `bowtie2 --version`;
	my ($bismark) = `bismark --version`;
	my ($bismark_genome_preparation) = `bismark_genome_preparation --version`;
	print "\n\n\e[1;33m ------------- BEGIN ------------ \e[0m\n";
	if (not defined $bedtools or $bedtools =~ /command not found/ or $bedtools =~ /bedtools v?([01].\d+|2\.0[0-9]|2\.1[0-6])/) {
		print "Please install bedtools at least version 2.17 before proceeding!\n";
		$bedtools = 0;
	}
	print "\n- \e[1;32m bedtools v2.17+ exists:\e[0m " . `which bedtools` if $bedtools ne 0;
	die if $bedtools eq 0;
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
	#	print "\n\n\e[1;33m ------------ BEGIN ------------> \e[0m\n";
}
use myFootLib; use FAlite;

#my $OPTS = "vp:"; getopts($OPTS);
#use vars   qw($opt_v $opt_p);
#my @VALS =   ($opt_v,$opt_p);
#my $MAINLOG = myFootLog::MAINLOG($0, \@VALS, $OPTS);

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

################
# ARGV Parsing #
###############

die "\nusage: $YW$0$N -n $CY<footPeak output folder>$N

${LGN}Options:$N
-G <gene>: Only run files with this gene in the name

" unless defined $opt_n and -d $opt_n;
my ($indexFile, $footPeakFolder) = ($opt_i, $opt_n);
my ($genewant) = $opt_G if defined $opt_G;

# sanity check -n footPeakFolder

($footPeakFolder) = getFullpath($footPeakFolder);
my $footClustFolder = "$footPeakFolder/FOOTCLUST/";
my $footKmerFolder  = "$footPeakFolder/KMER/";
my $uuid = getuuid();
my ($user) = $homedir =~ /home\/(\w+)/;
my $date = date();


############
# LOG FILE #
############
open (my $outLog, ">", "$footPeakFolder/footPeak_conv_logFile.txt") or die date() . ": Failed to create outLog file $footPeakFolder/footPeak_conv_logFile.txt: $!\n";
#open (my $outLog, ">", "$resDir/footPeak_conv_logFile.txt") or die "Failed to write to $resDir/footPeak_conv_logFile.txt: $!\n";
LOG($outLog, ">$0 version $version\n");
LOG($outLog, ">UUID: $uuid\n", "NA");
LOG($outLog, ">Date: $date\n", "NA");
LOG($outLog, ">Run script: $0 -n $opt_n\n", "NA");


##########
# OUTDIR #
##########
my $resDir = "$footPeakFolder";
	$resDir = getFullpath($resDir);

if (not -d $resDir) {
	makedir($resDir);
}

my ($OUTDIRS, $PEAK, $TEMP, $RCONV, $CPG, $ALL);
($OUTDIRS->{CONV}, $PEAK, $TEMP, $RCONV, $CPG, $ALL) = makeOutDir($resDir . "/CONVERSION/");

##############
# INDEX FILE #
##############


=comment
my %gene;
my @line = `cat $indexFile`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand, $feature) = split("\t", $line);
	$feature = "FEATURE_UNKNOWN" if not defined $feature;
	die "Undefine gene at line=$line\n" if not defined $gene;
	$gene{$gene}{feature} = $feature;
}
=cut
##########################
# PARSE FOOTPEAK LOGFILE #
##########################

my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
my $footLoop_run_script = `grep -iP "footLoop Run script\\s*:.+-g .+.fa" $footPeak_logFile`;
#LOG($outLog, "grep -iP \"footLoop Run script\\s*:.+-g .+.fa\" $footPeak_logFile\n");
DIELOG($outLog, "Cannot find footLoop_run_script $LCY$footLoop_run_script$N from footPeak logfile $footPeak_logFile\n") if not defined $footLoop_run_script or (defined $footLoop_run_script and $footLoop_run_script !~ /\w+/);
my @footLoop_run_script = split(" ", $footLoop_run_script);
my ($genomeFile, $xbuf);
for (my $i = 1; $i < @footLoop_run_script; $i++) {
	if ($footLoop_run_script[$i-1] eq "-g") {
		$genomeFile = $footLoop_run_script[$i];
	}
	if ($footLoop_run_script[$i-1] eq "-x") {
		$xbuf = $footLoop_run_script[$i];
	}
#	print "i=$i, $footLoop_run_script[$i-1]: $footLoop_run_script[$i]\n";
}
DIELOG($outLog, "Cannot find genome file from footPeak logfile $footPeak_logFile\n") if not defined $genomeFile or (defined $genomeFile and $genomeFile !~ /\w+/);
DIELOG($outLog, "Cannot find xbuf -x from footPeak logfile $footPeak_logFile\n") if not defined $xbuf or (defined $xbuf and $xbuf !~ /^\-?\d+$/);
LOG($outLog, "\ngenomeFile = $LCY$genomeFile$N\n");
LOG($outLog, "\nxbuf = $LGN$xbuf$N\n\n");

my %files;
my %coor;
LOG($outLog, "$footPeak_logFile\n");
my @lines   = `cat $footPeak_logFile`;
my ($label) = `cat $resDir/.LABEL`; chomp($label);
DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_logFile!\n\n") if not -e $footPeak_logFile;
DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $resDir/.LABEL!\n\n") if not -e "$resDir/.LABEL";
my ($thres, $window);
foreach my $line (@lines) {
	chomp($line);
	if ($line =~ /^[ \t]*def=.+, coor=.+/) {
		$line =~ s/(^\s+|\s+$)//g;
		my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
		if (defined $opt_G and $gene !~ /$opt_G/i) {
			LOG($outLog, date() . " Skipped $LCY$gene$N as it doesn't contain $LGN-G $opt_G$N\n");
			next;
		}
		$GENE = uc($GENE);
		DIELOG($outLog, "\n\ndied at processing $LCY$footPeak_logFile$N: can't parse index file def gene lqines\n\n$line\n\n") if not defined $STRAND;
		%{$coor{$GENE}} = ("CHR" => $CHR, "BEG" => $BEG, "END" => $END, "VAL" => $VAL);
		$coor{$GENE}{STRAND} = $STRAND eq "+" ? "Pos" : $STRAND eq "-" ? "Neg" : $STRAND =~ /^(Pos|Neg|Unk)$/ ? $STRAND : "Unk";
		$coor{$GENE}{BEG_BUFFER} = $BEG + $xbuf < 1 ? 1 : $BEG + $xbuf;
		LOG($outLog, "GENE=$GENE, STRAND=$STRAND, BEG=$BEG, xbuf=$xbuf\n");
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
my %notexists;
foreach my $GENE (sort keys %coor) {
	my $mygene = $GENE;
	my @strands = qw(Pos Neg Unk);
	my $geneStrand = $coor{$GENE}{STRAND};
	my ($CHR, $BEG_BUFFER) = ($coor{$GENE}{CHR}, $coor{$GENE}{BEG_BUFFER});
	LOG($outLog, "GENE=$GENE, CHR=$CHR, BEG+BUFFER = $BEG_BUFFER\n");
	#my $strand = $coor{$GENE}{STRAND};
	for (my $h1 = 0; $h1 < @strands; $h1++) {
		for (my $h2 = 0; $h2 < 4; $h2++) {
			my $readStrand = $strands[$h1];
			my $rconvType = $types[$h2];
			my $peakFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.out";
			my $nopkFile   = "$resDir/.CALL/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.out";
			$files{$peakFile} = $mygene;
			my ($flag) = getFlag($peakFile, $geneStrand, $readStrand, $rconvType);
			$flag =~ s/^PEAK//;
			LOG($outLog, "\n" . date() . ">GENE=$LRD$GENE$N $LGN$h1$LCY.$h2$N: $LPR$peakFile$N flag=$flag\n");
			my $outPEAKDir = $flag =~ /ALL/ ? $flag : "PEAK$flag";
			my $outNOPKDir = $flag =~ /ALL/ ? $flag : "NOPK$flag";
			LOG($outLog, date() . "$YW METHOD$LGN 1$N:\n");
			LOG($outLog, date() . " - $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_1peakonly.$outPEAKDir.tsv\n$N");
			LOG($outLog, date() . "$YW METHOD$LGN 2$N:\n");
			LOG($outLog, date() . " - $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_2allc.$outPEAKDir.tsv\n$N");
			LOG($outLog, date() . " - $LCY$resDir/CONVERSION/$outNOPKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.c_conv_2allc.$outNOPKDir.tsv\n$N");
			LOG($outLog, date() . "$YW METHOD$LGN 3$N:\n");
			LOG($outLog, date() . " - $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_3combined.$outPEAKDir.tsv\n$N");
			open (my $outPEAK_Method1peakonlyTSV, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_1peakonly.$outPEAKDir.tsv") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_1peakonly.$outPEAKDir.tsv$N: $!\n");
			open (my $outPEAK_Method2allcTSV, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_2allc.$outPEAKDir.tsv") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_2allc.$outPEAKDir.tsv$N: $!\n");
			open (my $outNOPK_Method2allcTSV, ">", "$resDir/CONVERSION/$outNOPKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.c_conv_2allc.$outNOPKDir.tsv") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outNOPKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.c_conv_2allc.$outNOPKDir.tsv$N: $!\n");
			open (my $outPEAK_Method3combinedTSV, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_3combined.$outPEAKDir.tsv") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_3combined.$outPEAKDir.tsv$N: $!\n");
			open (my $outPEAK_Method1peakonlyWIG, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_1peakonly.$outPEAKDir.wig") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_1peakonly.$outPEAKDir.wig$N: $!\n");
			open (my $outPEAK_Method2allcWIG, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_2allc.$outPEAKDir.wig") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_2allc.$outPEAKDir.wig$N: $!\n");
			open (my $outNOPK_Method2allcWIG, ">", "$resDir/CONVERSION/$outNOPKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.c_conv_2allc.$outNOPKDir.wig") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outNOPKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.NOPK.c_conv_2allc.$outNOPKDir.wig$N: $!\n");
			open (my $outPEAK_Method3combinedWIG, ">", "$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_3combined.$outPEAKDir.wig") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/CONVERSION/$outPEAKDir/$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType.PEAK.c_conv_3combined.$outPEAKDir.wig$N: $!\n");
			my @outpipeTSV = ($outPEAK_Method1peakonlyTSV, $outPEAK_Method2allcTSV, $outNOPK_Method2allcTSV, $outPEAK_Method3combinedTSV);
			my @outpipeWIG = ($outPEAK_Method1peakonlyWIG, $outPEAK_Method2allcWIG, $outNOPK_Method2allcWIG, $outPEAK_Method3combinedWIG);
			my @graphType  = qw(bar bar bar bar);
			my @methods = ("$outPEAKDir\t1peakonly", "$outPEAKDir\t2allc", "$outNOPKDir\t2allc", "$outPEAKDir\t3combined");
			my @colors = ("255,0,0", "0,155,0", "0,200,0", "155,0,155");
			my %res; my $max_x = -1;
			($res{PEAK}, $max_x) = calc_c_conv($peakFile, "PEAK", $rconvType, $flag, $resDir, $max_x, $outLog) if -e $peakFile and -s $peakFile > 0;
			($res{NOPK}, $max_x) = calc_c_conv($nopkFile, "NOPK", $rconvType, $flag, $resDir, $max_x, $outLog) if -e $nopkFile and -s $nopkFile > 0;
			if (not -e $peakFile or (-e $peakFile and -s $peakFile == 0)) {
				$notexists{$peakFile} = 1;
			}
			if (not -e $nopkFile or (-e $nopkFile and -s $nopkFile == 0)) {
				$notexists{$nopkFile} = 1;
			}
			my ($END_BUFFER) = $BEG_BUFFER + $max_x;
			my @XSPAN = ($BEG_BUFFER .. $END_BUFFER-1);
			for (my $i = 0; $i < @outpipeTSV; $i++) {
				print {$outpipeTSV[$i]} "label\tgene\treadStrand\twindow\tthres\trconvType\tpeakType\tmethods\t". join("\t", @XSPAN) . "\n";
				print {$outpipeTSV[$i]} "$label\t$mygene\t$readStrand\t$window\t$thres\t$rconvType\t$methods[$i]";
				print {$outpipeWIG[$i]} "track type=wiggle_0 graphType=$graphType[$i] name='$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType $methods[$i]' color=$colors[$i] visibility=2 maxHeightPixels=40:40:40 windowingFunction=mean viewLimits=0:1 autoScale=off\n";
				print {$outpipeWIG[$i]} "variableStep chrom=$CHR span=1\n";
			}

			my $printcount = 0;
			for (my $i = 0; $i < $max_x; $i++) {
				my $pos = $BEG_BUFFER + $i;
				LOG($outLog, " --> i = $i, BEG_BUFFER = $BEG_BUFFER, pos = beg-buffer+i = $pos\n") if $i % 100 == 0;
				my $peakconv1 = defined $res{PEAK}->[$i]{1}{total} ? int(1000 * $res{PEAK}->[$i]{1}{total} / $res{PEAK}->[$i]{1}{count} + 0.5)/1000 : "NA";

				print "\tres{NOPK} defined but undef at peakfile=$LCY$peakFile undef res{PEAK} i=$i 2\n" if $printcount == 0 and not defined $res{PEAK}->[$i]{2} and defined $res{NOPK}->[$i]{2};
				print "\tres{PEAK} defined but undef at nopkfile=$LCY$nopkFile undef res{NOPK} i=$i 2\n" if $printcount == 0 and not defined $res{NOPK}->[$i]{2} and defined $res{PEAK}->[$i]{2};
				LOG($outLog, date() . "\tres{NOPK} defined but undef at peakfile=$LCY$peakFile undef res{PEAK} i=$i 2\n", "NA") if not defined $res{PEAK}->[$i]{2} and defined $res{NOPK}->[$i]{2};
				LOG($outLog, date() . "\tres{PEAK} defined but undef at nopkfile=$LCY$nopkFile undef res{NOPK} i=$i 2\n", "NA") if not defined $res{NOPK}->[$i]{2} and defined $res{PEAK}->[$i]{2};
				my $peakconv2 = defined $res{PEAK}->[$i]{2}{total} ? int(1000 * $res{PEAK}->[$i]{2}{total} / $res{PEAK}->[$i]{2}{count} + 0.5)/1000 : "NA";
				my $nopkconv2 = defined $res{NOPK}->[$i]{2}{total} ? int(1000 * $res{NOPK}->[$i]{2}{total} / $res{NOPK}->[$i]{2}{count} + 0.5)/1000 : "NA";

				my $peakconv3 = ($peakconv2 eq "NA" and $nopkconv2 eq "NA") ? "NA" :
									 ($peakconv2 eq "NA") ? $nopkconv2 :
									 ($nopkconv2 eq "NA") ? $peakconv2 :
									 int(1000 * ($res{PEAK}->[$i]{2}{total} + $res{NOPK}->[$i]{2}{total} ) / ($res{PEAK}->[$i]{2}{count} + $res{NOPK}->[$i]{2}{count} ) + 0.5)/1000;
				print $outPEAK_Method1peakonlyTSV "\t$peakconv1";
				print $outPEAK_Method2allcTSV "\t$peakconv2";
				print $outNOPK_Method2allcTSV "\t$nopkconv2";
				print $outPEAK_Method3combinedTSV "\t$peakconv3";
				print $outPEAK_Method1peakonlyWIG "$pos\t$peakconv1\n" if $peakconv1 ne "NA";
				print $outPEAK_Method2allcWIG "$pos\t$peakconv2\n" if $peakconv2 ne "NA";
				print $outNOPK_Method2allcWIG "$pos\t$nopkconv2\n" if $nopkconv2 ne "NA";
				print $outPEAK_Method3combinedWIG "$pos\t$peakconv3\n" if $peakconv3 ne "NA";
				$printcount = 1;
			}
			for (my $i = 0; $i < @outpipeTSV; $i++) {
				print {$outpipeTSV[$i]} "\n";
			}
		}
	}
}

LOG($outLog, "\n\n\n$YW -------------- FILES THAT DIDN'T EXISTS -------------- $N \n\n\n");
foreach my $notexist (sort keys %notexists) {
	LOG($outLog, " - $notexist\n");
}
###############
# SUBROUTINES #
###############

# 0 is bad or not data (6)
# 1 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak

sub calc_c_conv {
	my ($file, $peaktype, $rconvType, $flag, $resDir, $max_x, $outLog) = @_;
	my $isCG = 0;
	if ($rconvType =~ /(CG|GC)/ or $flag =~ /_C/) {
		$isCG = 1;
	}
	my @res;
	open (my $in, "<", $file) or DIELOG($outLog, date() . " FailedD to read from $file: $!\n");
	while (my $line = <$in>) {
		chomp($line);
		my ($read, @conv) = split("\t", $line);
		for (my $i = 0; $i < @conv; $i++) {
			next if $conv[$i] < 4; # 0 or 1
			next if $isCG eq 0 and $conv[$i] =~ /^(5|7|9)$/;
			if ($peaktype =~ /PEAK/) {
				$res[$i]{1}{total} += $conv[$i] >= 8 ? 1 : 0; # method 1 (peak only)
				$res[$i]{1}{count} ++;
			}
			$res[$i]{2}{total} += $conv[$i] >= 6 ? 1 : 0; # method 2 (all c regardless of peak)
			$res[$i]{2}{count} ++;
			$max_x = $i if $max_x < $i;
		}
	}
	close $in;
	return (\@res, $max_x);
}

__END__
###############################
# Get peaks from PEAKS_GENOME #
###############################
my $cluster;
my @bedFiles = <$footPeakFolder/PEAKS_GENOME/*.genome.bed>;
my @files;
LOG($outLog, date() . "1. Getting cluster and preprocessing bed files!\n");
foreach my $bedFile (sort @bedFiles) {
	if (defined $opt_G) {next if $bedFile !~ /$opt_G/};
	my ($bedFileName) = $bedFile =~ /^$footPeakFolder\/PEAKS_GENOME\/(.+.PEAK.genome.bed)$/;
	my $tempFile;
	my @alphabets = qw(A B C W D E F);
	my @treats = qw(dens cont skew);
	for (my $i = 0; $i < @alphabets; $i++) {
		for (my $j = 0; $j < @treats; $j++) {
			$tempFile = "$footPeakFolder/GCPROFILE/$bedFileName";
			my $tsvFile = "$footPeakFolder/GCPROFILE/$bedFileName\_100\_$alphabets[$i].temp.fa.$treats[$j].tsv";
			#print "\t- $LCY$tsvFile$N\n";
			@files = (@files, $tsvFile);
		}
	}
	$tempFile =~ s/\/+/\//g;
	$cluster = get_cluster($bedFile, $tempFile, $cluster, $outLog);
	#print "$LCY$bedFile$N\n";
	preprocess_bed($bedFile, $outLog) if not defined $opt_S;
}
LOG($outLog, date() . "2. Calculating GC skew (might take a couple minutes)\n");
if (not defined ($opt_S)) {
	system("run_script_in_paralel2.pl \"fastaFromBed -fi $genomeFile -bed FILENAME -fo FILENAME.fa -s -name\" $resDir.TEMP/ \"_[ABCDEFW].temp\" 1 > $resDir/.TEMP/fastaFromBed.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run fastaFromBed: $!\n");
	system("run_script_in_paralel2.pl \"rename.pl FILENAME PCB .PCB\" $resDir.TEMP/ temp 1  > $resDir/.TEMP/rename.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run rename.pl: $!\n");
	system("run_script_in_paralel2.pl \"counter_cpg_indiv.pl -w 200 -s 1 -o $resDir -A FILENAME\" $resDir.TEMP/ _100.+temp.fa 1  > $resDir/.TEMP/counter_cpg_indiv.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run counter_cpg_indiv.pl: $!\n");
}
####### PARAMETERS
sub get_cluster {
	my ($bedFile, $tempFile, $cluster, $outLog) = @_; 
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my ($clusterFile) = "$footPeakFolder/FOOTCLUST/.TEMP/$bedFilename";
	$clusterFile =~ s/.genome.bed/.local.bed.clust/;
	my $maxClust = -1;
	if (-e $clusterFile) {
		my $linecount = -1;
		LOG($outLog, date() . "   bedFile=$LCY$bedFile$N, clusterFile=$LGN$clusterFile$N\n","NA");
		open (my $clusterIn, "<", $clusterFile) or DIELOG($outLog, "Failed to read from clusterFile $clusterFile: $!\n");
		while (my $line = <$clusterIn>) {
			chomp($line);
			$linecount ++;
			next if $linecount == 0;
			my ($id, $x, $xmax, $y, $ymax, $clust) = split("\t", $line);
			$maxClust = $clust if $maxClust < $clust;
			my ($id2, $number) = $id =~ /^(\d+)\.(\d+)$/;
			($id2) = $id if not defined $id2;
#			print "id=$id, id2=$id2, number=$number, $line\n" if $linecount < 10;
			$id = $id2;
			if (not defined $cluster->{$tempFile}{$id}) {
				my $currlen = $xmax - $x;
				#print "$clusterFile id=$id clust=$clust num=$number len=$currlen\n$tempFile,$id,$clust\n" if $id eq "1704132300151533640";
				$cluster->{$tempFile}{$id}{clust} = $clust;
				$cluster->{$tempFile}{$id}{len} = ($xmax - $x);
			}
			elsif (defined $cluster->{$tempFile}{$id} and $cluster->{$tempFile}{$id}{len} < $xmax - $x) {
				my $currlen = $xmax - $x;
				LOG($outLog, date() . " $clusterFile id=$id num=$number len=$cluster->{$tempFile}{$id}{len} < $currlen\n","NA");# if $id eq "1704132300151533640";
				$cluster->{$tempFile}{$id}{clust} = $clust;
				$cluster->{$tempFile}{$id}{len} = ($xmax - $x);
			}
		}
	}
	LOG($outLog, date() . "$LCY$bedFile$N cluster=$LGN$maxClust$N\n","NA");
	#print "temp=$LPR$tempFile$N\n";
	return($cluster);
}

sub preprocess_bed {
	my ($bedFile, $outLog) = @_; 
	LOG($outLog, date() . " Getting fasta and calculating GC skew from: $LCY$bedFile$N\n", "NA");
	my ($bedFolder, $bedFilename) = getFilename($bedFile, "folderfull");
	my $window = 100;
	my $window2 = 200;
	my $outputA = "$resDir/.TEMP/$bedFilename\_$window\_A.temp";
	my $outputB = "$resDir/.TEMP/$bedFilename\_$window\_B.temp";
	my $outputC = "$resDir/.TEMP/$bedFilename\_$window\_C.temp";
	my $outputD = "$resDir/.TEMP/$bedFilename\_$window\_D.temp";
	my $outputE = "$resDir/.TEMP/$bedFilename\_$window\_E.temp";
	my $outputF = "$resDir/.TEMP/$bedFilename\_$window\_F.temp";
	my $outputW = "$resDir/.TEMP/$bedFilename\_$window\_W.temp";

	system("bedtools_bed_change.pl -a -x -$window2 -y 0 -i $bedFile -o $outputA > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x -$window -y $window -i $bedFile -o $outputB > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -a -x 0 -y $window2 -i $bedFile -o $outputC > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -x 0 -y 0 -i $bedFile -o $outputW > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window2 -y 0 -i $bedFile -o $outputD > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x -$window -y $window -i $bedFile -o $outputE > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
	system("bedtools_bed_change.pl -b -x 0 -y $window2 -i $bedFile -o $outputF > $outputA.LOG 2>&1") == 0 or DIELOG($outLog, "Failed to run bedtools_bed_change.pl: $!\n");
}

#my @files;# = <$resDir/PCB*.tsv>;
#if (defined $opt_G) {
#	@files = <$resDir/PCB*$opt_G*.tsv>;
#}
#else {
#	@files = <$resDir/PCB*.tsv>;
#}
LOG($outLog, date() . "3. Processing " . scalar(@files) . " files in $LCY$resDir$N\n");
my %data;
my @header = ("label", "gene", "strand", "window", "threshold", "convtype", "wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$resDir/$N!\n" and exit if (@files == 0);
foreach my $input1 (sort @files) {
	my ($tempFile) = $input1 =~ /^(.+)_100_.\.temp.fa.\w+.tsv$/;
	$tempFile =~ s/\/+/\//g;
	#print "input1=$input1, tempFile=$LCY$tempFile$N\n";
	$input1 =~ s/\/+/\//g;
	my ($WINDOW, $SAMPLE, $TYPE);
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	if (defined $opt_G) {
		next if $fileName1 !~ /$opt_G/;
	}
#	print "$input1\n";
	#my ($label, $barcode, $desc, $gene, $strand, $window, $threshold, $convtype, $wind2, $sample, $type) = $fileName1 =~ /^(PCB.+)_(BC\d+)?_?(\w+)?_?gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	next if $fileName1 !~ /^PCB/;
	my @arr = $fileName1 =~ /^(PCB.+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	$arr[0] =~ s/^(PCB\d+)_.+$/$1/;
	if (not defined $arr[0]) {
		for (my $i = 0; $i < @arr; $i++) {
			DIELOG($outLog, date() . "fileName=$fileName1 Undefined i=$i header=$header[$i] arr[i] undef\n") if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		}
	}
	my $outName = join("_", @arr[0..6]) . "_" . $arr[8];
	for (my $i = 0; $i < @arr; $i++) {
		DIELOG($outLog, date() . "Undefined i=$i header=$header[$i] arr[i] undef\n") if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		$data{data}{$outName}{$header[$i]} = $arr[$i];
		$WINDOW = $arr[$i] if $header[$i] eq "wind2";
		$SAMPLE = $arr[$i] if $header[$i] eq "sample";
		$TYPE = $arr[$i] if $header[$i] eq "type";
	}
	LOG($outLog, date() . "  - gene=$LGN$arr[1]$N pos=$LGN$arr[7]$N type=$LPR$arr[8]$N input=$outName$N\n","NA") if $arr[8] eq "skew";
	open (my $in1, "<", $input1) or DIELOG($outLog, date() . " Cannot read from $input1: $!\n");
	my $linecount = 0;
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		$linecount ++;
		my ($read, $value) = split("\t", $line);
		my ($id1, $id2, $id3) = $read =~ /^.*m(\d+_\d+)_\d+_c\w+_\w+_\w+\/(\d+)\/(ccs|\d+_\d+)/;
		$id3 = 0 if not defined $id3;
		$id3 = 0 if $id3 eq "ccs";
		DIELOG($outLog, "Failed to parse id1/2/3 from read=$LCY$read$N, file=$LGN$input1$N\n") if not defined $id1 or not defined $id2 or not defined $id3;
		my $id = "$id1$id2$id3"; $id =~ s/_//g;
		my $cluster = $cluster->{$tempFile}{$id}{clust}; $cluster = -1 if not defined $cluster;
		#print "id=$id, cluster=$cluster\n$input1,$id,$cluster\n" if $id eq "1704132300151533640";
		$data{cluster}{$outName}{$read} = $cluster;
		$data{id}{$outName}{$read} = $id;
		$data{read}{$outName}{$read}{$SAMPLE} = $value;
		$data{input}{$outName}{$SAMPLE} = join("_", @arr[0..5]);
	}	
	close $in1;
}
open (my $out1, ">", "$resDir/RESULT.TSV") or DIELOG($outLog, date() . "Cannot write to $resDir/RESULT.TSV: $!\n");
foreach my $outName (sort keys %{$data{input}}) {
	#open (my $out1, ">", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	print $out1 "file\tid\tcluster\tread\twindow\ttype\tfeature";
	foreach my $sample (sort keys %{$data{input}{$outName}}) {
		print $out1 "\t$sample";
	}
	print $out1 "\n";
	last;
	#close $out1;
}

foreach my $outName (sort keys %{$data{read}}) {
#open (my $out1, ">>", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	my $GENE = $data{data}{$outName}{gene};
	my $feature = $gene{$GENE}{feature}; 
	$feature = "FEATURE_UNKNOWN" if not defined $feature;#print "Undef gene=$GENE feature\n" and next if not defined $feature;
	foreach my $read (sort keys %{$data{read}{$outName}}) {
		print $out1 "$data{input}{$outName}{$SAMPLE}\t$data{id}{$outName}{$read}\t$data{cluster}{$outName}{$read}\t$read\t$WINDOW\t$TYPE\t$feature";
		foreach my $sample (sort keys %{$data{read}{$outName}{$read}}) {
			print $out1 "\t$data{read}{$outName}{$read}{$sample}";
		}
		print $out1 "\n";
	}
	#close $out1;
}

GCprofile_Rscript($resDir, $outLog);

sub GCprofile_Rscript {
	my ($resDir, $outLog) = @_;
my $resDirFullpath = getFullpath($resDir);

my $RESULT = $resDir . "/RESULT.TSV";
my $LABEL = `cat $resDir/../.LABEL`; DIELOG($outLog, "Cannot find $resDir/../.LABEL!\n") if not defined $LABEL; chomp($LABEL);
$LABEL = $resDirFullpath . "/$LABEL";
my $Rscript = "
.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\",
\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\",
.libPaths()) )
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)    
library(ggplot2);
library(reshape2);

RESULT=\"$RESULT\";
PDFCLUSTCpGdens = \"$LABEL\_BYCLUST_CpGdens.pdf\"
PDFGENESCpGdens = \"$LABEL\_BYGENES_CpGdens.pdf\"
PDFCLUSTGCcont = \"$LABEL\_BYCLUST_GCcont.pdf\"
PDFGENESGCcont = \"$LABEL\_BYGENES_GCcont.pdf\"
PDFCLUSTGCskew = \"$LABEL\_BYCLUST_GCskew.pdf\"
PDFGENESGCskew = \"$LABEL\_BYGENES_GCskew.pdf\"
PDFCLUST = c(PDFCLUSTCpGdens,PDFCLUSTGCcont,PDFCLUSTGCskew)
PDFGENES = c(PDFGENESCpGdens,PDFGENESGCcont,PDFGENESGCskew)

df = read.table(RESULT,header=T,sep=\"\\t\")
dm = melt(df,id.vars=c(\"file\",\"id\",\"cluster\",\"read\",\"window\",\"type\",\"feature\"));
dm\$feature = factor(dm\$feature,levels=c(\"PROMOTER\",\"GENEBODY\",\"TERMINAL\",\"FEATURE_UNKNOWN\"));
dm\$variable = factor(dm\$variable,levels=c(\"A\",\"B\",\"C\",\"W\",\"D\",\"E\",\"F\"))
dm\$gene = paste(dm\$file,dm\$feature)
genes = unique(dm\$gene)
dm\$group = paste(dm\$cluster)
##
clustTot = dim(aggregate(dm\$window,by=list(dm\$gene,dm\$cluster),sum))[1]
genesTot = dim(aggregate(dm\$window,by=list(dm\$gene),sum))[1]
clustCount = as.data.frame(plyr::count(dm,c(\"cluster\",\"gene\")));
clustCount\$freq = clustCount\$freq / (length(unique(dm\$variable)) * length(unique(dm\$type))); colnames(clustCount)[3] = \"clustGroup\"
genesCount = as.data.frame(aggregate(clustCount\$clustGroup,by=list(clustCount\$gene),sum));colnames(genesCount) = c(\"gene\",\"genesGroup\");
dm = merge(dm,clustCount,by=c(\"cluster\",\"gene\"),all=T)
dm = merge(dm,genesCount,by=c(\"gene\"),all=T)
dm\$clustGroup = paste(dm\$file,\" (\",dm\$feature,\") cluster \",dm\$cluster,\" (\",dm\$clustGroup,\" reads)\",sep=\"\")
dm\$genesGroup = paste(dm\$file,\" (\",dm\$feature,\") (\",dm\$genesGroup,\" reads)\",sep=\"\")
##
types = c(\"dens\",\"cont\",\"skew\")
ylimsMin=c(0,0,-1)
ylimsMax=c(1.2,1,1)
ylines=c(0.6,0.5,0)
ylabs = c(\"CpG Density\",\"GC Content\",\"GC Skew\")
for (i in 1:length(PDFCLUST)) {
pdf(PDFCLUST[i],width=7,height=7*clustTot)
temp = dm[dm\$type == types[i],]
p = ggplot(temp,aes(variable,value)) +
geom_boxplot(aes(fill=variable),outlier.shape=NA) +
theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin[i],ylimsMax[i])) +
annotate(geom=\"segment\",x=0,xend=8,y=ylines[i],yend=ylines[i],lty=2) +
facet_grid(clustGroup~.) +
ylab(ylabs[i]) + xlab(\"Samples\")
print(p)
dev.off()
}

for (i in 1:length(PDFGENES)) {
pdf(PDFGENES[i],width=7,height=7*genesTot)
temp = dm[dm\$type == types[i],]
p = ggplot(temp,aes(variable,value)) +
geom_boxplot(aes(fill=variable),outlier.shape=NA) +
theme_bw() + theme(panel.grid=element_blank(),legend.position=\"none\") + coord_cartesian(ylim=c(ylimsMin[i],ylimsMax[i])) +
annotate(geom=\"segment\",x=0,xend=8,y=ylines[i],yend=ylines[i],lty=2) +
facet_grid(genesGroup~.) +
ylab(ylabs[i]) + xlab(\"Samples\")
print(p)
dev.off()
}
";

LOG($outLog, date() . "4. $YW Running R script $LCY$resDir/RESULT.R$N\n");
open (my $outR, ">", "$resDir/RESULT.R") or DIELOG($outLog, date() . " Failed to write to $LCY$resDir/RESULT.R$N: $!\n");
print $outR $Rscript;
close $outR;
system("run_Rscript.pl $resDir/RESULT.R > $resDir/RESULT.R.LOG 2>&1") == 0 or DIELOG($outLog, date() . " Failed to run $LCY$resDir/RESULT.R$N: $!\n");
my $tailR = `tail $resDir/RESULT.R.LOG`;

LOG($outLog, "\n\n" . date() . " ${LGN}SUCCESS on running $LCY$resDir/RESULT.R$YW.\nLast 5 rows of log message:$N\n$N$tailR


${YW}To Run R script manually, do:$N
run_Rscript.pl $resDir/RESULT.R


${YW}Outputs:$N

- $LGN#grouped by each gene:$N
$LCY$LABEL\_BYGENES_CpGdens.pdf$N
$LCY$LABEL\_BYGENES_GCcont.pdf$N
$LCY$LABEL\_BYGENES_GCskew.pdf$N

- $LGN#grouped by each gene and each cluster:$N
$LCY$LABEL\_BYCLUST_CpGdens.pdf$N
$LCY$LABEL\_BYCLUST_GCcont.pdf$N
$LCY$LABEL\_BYCLUST_GCskew.pdf$N


");
}

__END__



close $out1;

__END__
PCB1_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed_100_E.temp.fa.dens.tsv
__END__
#				print $outpipeTSV[$i] "$label\_gene$mygene\_$readStrand\_$window\_$thres\_$rconvType\t$methods[$i]";

