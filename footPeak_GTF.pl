#!/usr/bin/perl

use warnings;use strict;use Getopt::Std; use Cwd qw(abs_path);use File::Basename qw(dirname);
use vars qw($opt_i $opt_f $opt_x $opt_y $opt_o $opt_p $opt_t $opt_n $opt_G $opt_v);
getopts("i:f:x:y:o:p:t:n:G:v");

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
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

my ($footPeakFolder) = ($opt_n);

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$LGN [Optional: -G <gene to process>]$N -n $LCY<footPeak output folder>$N

";

die $usage if not defined $opt_n;

my $outFolder = "$footPeakFolder/GTF";
makedir($outFolder);

my ($footPeak_logFile) = "$footPeakFolder/footPeak_logFile.txt";
my ($logFile) = "$footPeakFolder/footPeakGTF_logFile.txt";
open (my $outLog, ">", $logFile) or print "\n\nFailed to write to $logFile: $!\n\n" and die;
print "LOGFILE:\n$YW$logFile$N\n";
my $opts = parse_GTF_footPeak_logFile($footPeakFolder, $footPeak_logFile, $outLog);

my ($bufferL) = $opts->{footLoop2}{x};
my ($bufferR) = $opts->{footLoop2}{y};
die "Left (-x $bufferL) and Right (-y $bufferR) buffer must be integer!\n" if $bufferL !~ /^\-?\d+$/ or $bufferR !~ /^\-?\d+$/;
#$opt_p
#$opt_o
#$opt_t
#$opt_i
#$opt_f
my ($footLoopFolder) = $opts->{footLoop2}{n}; 
DIELOG($outLog, "UNDEFINED footLoopFolder\n") if not defined $footLoopFolder;
DIELOG($outLog, "Does not exist footLoopFolder=$LCY$footLoopFolder$N\n") if not -e $footLoopFolder;
my $coor = parse_geneIndexFile($footLoopFolder, $outLog);

my (@callFiles) = <$footPeakFolder/.CALL/*.out>;

my $label = "";
if (-e "$footPeakFolder/.LABEL") {
	print "cat $footPeakFolder/.LABEL";
   ($label) = `cat $footPeakFolder/.LABEL`;
   chomp($label);
}
else {
   DIELOG($outLog, "Failed to parse label from .LABEL in $footPeakFolder/.LABEL\n");
}

foreach my $callFile (sort @callFiles) {
   if (defined $opt_G and $callFile !~ /$opt_G/i) {
      LOG($outLog, date() . " Skipped $LCY$callFile$N as it doesn't contain $LGN-G $opt_G$N\n");
      next;
   }

	my ($folder1, $fileName1) = getFilename($callFile, "folderfull");
   # get gene and strand from file name
   my $parseName = parseName($fileName1);
   my ($label2, $gene, $strand, $window, $thres, $type) = @{$parseName->{array}};
   LOG($outLog, "Using label=$label2. Inconsistent label in filename $LCY$fileName1$N\nLabel from $footPeakFolder/.LABEL: $label\nBut from fileName: $label2\n\n") if $label2 ne $label;
   $label = $label2;
#	next unless $label eq "PCB7" and $strand eq "Pos" and $gene eq "RPS24" and $type eq "CH";
#	my ($label, $gene, $strand, $window, $thres, $type) = $fileName1 =~ /^(.+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(GC|GH|CG|CH)/;
	my $STRAND = $coor->{$gene}{strand};
	my $clustFile = $callFile =~ /.PEAK.out/ ? "$footPeakFolder/FOOTCLUST/CLUST_LOCAL/$label\_gene$gene\_$strand\_$window\_$thres\_$type.PEAK.local.bed.indiv.clust" : "";

#CLUST_LOCAL/PCB13_bcBC3_plasmidPFC8_descSUPERCOILED_genePFC8_SNRPN_REVERSE_Neg_20_0.55_GC.PEAK.local.bed.indiv.clust
#.TEMP/PCB13_bcBC3_plasmidPFC8_descSUPERCOILED_genePFC8_SNRPN_REVERSE_Neg_20_0.55_GC.PEAK.local.bed.clust

	my @folderz = qw(PEAK PEAK_C PEAK_TEMP PEAK_TEMP_C PEAK_RCONV PEAK_RCONV_C PEAK_TEMP_RCONV PEAK_TEMP_RCONV_C NOPK NOPK_C NOPK_TEMP NOPK_TEMP_C NOPK_RCONV NOPK_RCONV_C NOPK_TEMP_RCONV NOPK_TEMP_RCONV_C ALL);
	foreach my $folderz (@folderz) {
		makedir("$outFolder/$folderz") if not -d "$outFolder/$folderz";
	}
	my $fwstrand = $STRAND eq "+" ? "Pos" : "Neg";
	my $rvstrand = $STRAND eq "+" ? "Neg" : "Pos";
	my $fwconv 	 = $STRAND eq "+" ? "CH" : "GH";
	my $rvconv   = $STRAND eq "+" ? "GH" : "CH";
	my $fwconv_c = $STRAND eq "+" ? "CG" : "GC";
	my $rvconv_c = $STRAND eq "+" ? "GC" : "CG";
	my $gtfFolder = "ALL";
	$gtfFolder = "PEAK" if $callFile =~ /$fwstrand.+$fwconv.PEAK.out/;
	$gtfFolder = "NOPK" if $callFile =~ /$fwstrand.+$fwconv.NOPK.out/;
	$gtfFolder = "PEAK_C" if $callFile =~ /$fwstrand.+$fwconv_c.PEAK.out/;
	$gtfFolder = "NOPK_C" if $callFile =~ /$fwstrand.+$fwconv_c.NOPK.out/;
	$gtfFolder = "PEAK_TEMP" if $callFile =~ /$rvstrand.+$rvconv.PEAK.out/;
	$gtfFolder = "NOPK_TEMP" if $callFile =~ /$rvstrand.+$rvconv.NOPK.out/;
	$gtfFolder = "PEAK_TEMP_C" if $callFile =~ /$rvstrand.+$rvconv_c.PEAK.out/;
	$gtfFolder = "NOPK_TEMP_C" if $callFile =~ /$rvstrand.+$rvconv_c.NOPK.out/;
	$gtfFolder = "PEAK_RCONV" if $callFile =~ /$fwstrand.+$rvconv.PEAK.out/;
	$gtfFolder = "NOPK_RCONV" if $callFile =~ /$fwstrand.+$rvconv.NOPK.out/;
	$gtfFolder = "PEAK_RCONV_C" if $callFile =~ /$fwstrand.+$rvconv_c.PEAK.out/;
	$gtfFolder = "NOPK_RCONV_C" if $callFile =~ /$fwstrand.+$rvconv_c.NOPK.out/;
	$gtfFolder = "PEAK_TEMP_RCONV" if $callFile =~ /$rvstrand.+$fwconv.PEAK.out/;
	$gtfFolder = "NOPK_TEMP_RCONV" if $callFile =~ /$rvstrand.+$fwconv.NOPK.out/;
	$gtfFolder = "PEAK_TEMP_RCONV_C" if $callFile =~ /$rvstrand.+$fwconv_c.PEAK.out/;
	$gtfFolder = "NOPK_TEMP_RCONV_C" if $callFile =~ /$rvstrand.+$fwconv_c.NOPK.out/;

#	my ($peak) = $callFile =~ /$fwstrand.+$fwconv.PEAK.out/ ? "PEAK" : $callFile =~ /$fwstrand.+$fwconv_c.PEAK.out/ ? "PEAK_C" :
#					 $callFile =~ /$fwstrand.+CG.NOPK.out/ ? "NOPK_C" : $callFile =~ /$fwstrand.+CH.NOPK.out/ ? "NOPK" :
#					 $callFile =~ /$fwstrand.+GC.PEAK.out/ ? "PEAK_TEMP_C" : $callFile =~ /$fwstrand.+GH.PEAK.out/ ? "PEAK_TEMP" :
##					 $callFile =~ /$fwstrand.+GC.NOPK.out/ ? "NOPK_TEMP_C" : $callFile =~ /$fwstrand.+GH.NOPK.out/ ? "NOPK_TEMP" :
#					 $callFile =~ /$fwstrand.+CG.PEAK.out/ ? "PEAK_C" : $callFile =~ /$fwstrand.+CH.PEAK.out/ ? "PEAK" :
#					 $callFile =~ /$fwstrand.+CG.NOPK.out/ ? "NOPK_C" : $callFile =~ /$fwstrand.+CH.NOPK.out/ ? "NOPK" :
#					 $callFile =~ /$fwstrand.+GC.PEAK.out/ ? "PEAK_TEMP_C" : $callFile =~ /$fwstrand.+GH.PEAK.out/ ? "PEAK_TEMP" :
#					 $callFile =~ /$fwstrand.+GC.NOPK.out/ ? "NOPK_TEMP_C" : $callFile =~ /$fwstrand.+GH.NOPK.out/ ? "NOPK_TEMP" :
#
	my $outName = "$label\_gene$gene\_$strand\_$window\_$thres\_$type.$gtfFolder";
	my $output = "$outFolder/$gtfFolder/$outName.gtf";
	my $outputColor = "$LCY$outFolder/LGN$gtfFolder/$LPR$outName.gtf$N";
#	if ($peak eq "PEAK") {
#		$output = "$outFolder/PEAK/$outName.gtf" if $strand eq "Pos" and $STRAND eq "+" and $type eq "CG" and $gene =~ /FORWARD/;
#		$output = "$outFolder/PEAK/$outName.gtf" if $strand eq "Neg" and $STRAND eq "-" and $type eq "GC" and $gene =~ /REVERSE/;
#		$output = "$outFolder/PEAK/$outName.gtf" if $strand eq "Pos" and $STRAND eq "+" and $type eq "CH" and $gene !~ /FORWARD/;
#		$output = "$outFolder/PEAK/$outName.gtf" if $strand eq "Neg" and $STRAND eq "-" and $type eq "GH" and $gene !~ /REVERSE/;
#		$output = "$outFolder/PEAKNEG/$outName.gtf" if $strand eq "Pos" and $STRAND eq "-" and $type eq "CG" and $gene =~ /REVERSE/;
#		$output = "$outFolder/PEAKNEG/$outName.gtf" if $strand eq "Neg" and $STRAND eq "+" and $type eq "GC" and $gene =~ /FORWARD/;
#		$output = "$outFolder/PEAKNEG/$outName.gtf" if $strand eq "Pos" and $STRAND eq "-" and $type eq "CH" and $gene !~ /REVERSE/;
#		$output = "$outFolder/PEAKNEG/$outName.gtf" if $strand eq "Neg" and $STRAND eq "+" and $type eq "GH" and $gene !~ /FORWARD/;
#	}
#	elsif ($peak eq "NOPK") {
#		$output = "$outFolder/NOPK/$outName.gtf" if $strand eq "Pos" and $STRAND eq "+" and $type eq "CG" and $gene =~ /FORWARD/;
#		$output = "$outFolder/NOPK/$outName.gtf" if $strand eq "Neg" and $STRAND eq "-" and $type eq "GC" and $gene =~ /REVERSE/;
#		$output = "$outFolder/NOPK/$outName.gtf" if $strand eq "Pos" and $STRAND eq "+" and $type eq "CH" and $gene !~ /FORWARD/;
#		$output = "$outFolder/NOPK/$outName.gtf" if $strand eq "Neg" and $STRAND eq "-" and $type eq "GH" and $gene !~ /REVERSE/;
#		$output = "$outFolder/NOPKNEG/$outName.gtf" if $strand eq "Pos" and $STRAND eq "-" and $type eq "CG" and $gene =~ /REVERSE/;
#		$output = "$outFolder/NOPKNEG/$outName.gtf" if $strand eq "Neg" and $STRAND eq "+" and $type eq "GC" and $gene =~ /FORWARD/;
#		$output = "$outFolder/NOPKNEG/$outName.gtf" if $strand eq "Pos" and $STRAND eq "-" and $type eq "CH" and $gene !~ /REVERSE/;
#		$output = "$outFolder/NOPKNEG/$outName.gtf" if $strand eq "Neg" and $STRAND eq "+" and $type eq "GH" and $gene !~ /FORWARD/;
#	}
#	next if $output =~ /ALL/;
	LOG($outLog, "\n\nOutput=$outputColor\n");
	my $color = $strand eq "Pos" ? "200,50,0" : $strand eq "Neg" ? "0,50,200" : "50,155,50";
	my $header = "track name=$outName color=$color";
	turn_into_gtf($gene, $callFile, $coor->{$gene}, $output, $header, $outLog, $clustFile, $gtfFolder) if -s $callFile != 0;
#	parse($gene, $negFile, $data{$gene}, $GTFNegfile) if -e $negFile and -s $negFile != 0;
#	print "\tgene $gene\'s neg file $negFile does not exist or has nothign in it!\n" if not -e $negFile or (-e $negFile and -s $negFile == 0);
}

close $outLog;

sub turn_into_gtf {
	my ($gene, $callFile, $coor, $output, $header, $outLog, $clustFile, $gtfFolder) = @_;
	my %clust;
	my %clustlen;
	if ($clustFile ne "" and -e $clustFile) {
		open (my $clustFileIn, "<", $clustFile) or DIELOG($outLog, "Failed to read from $clustFile: $!\n");
		my $linecount = 0;
		my ($total_line) = `wc -l $clustFile` =~ /^\s*(\d+)/;
		while (my $line = <$clustFileIn>) {
			chomp($line); $linecount ++;
			next if $line =~ /^id/;
			my ($id2, $x, $xmax, $y, $ymax, $clust, $clust2) = split("\t", $line);
			my $id = $id2;
			my $idcount = 0;
			($id, $idcount) = $id =~ /^(.+)(\d)$/;
			$id = $id2;
			my $idlen = $xmax - $x;
			if (not defined $clustlen{$id} or (defined $clustlen{$id} and $clustlen{$id} < $idlen)) {
				$clust{$id} = $total_line > 500 ? int($y/$total_line * 500) : $y;
				$clustlen{$id} = $idlen;
			}
		}
		print "$LGN exist$N $clustFile!\n";
	}
	else {
		print "$LRD not exist$N $clustFile\n";
	}
	my ($CHR, $BEG, $END, $STRAND) = ($coor->{chr}, $coor->{beg}, $coor->{end}, $coor->{strand});
	$STRAND = $STRAND eq "Pos" ? "+" : $STRAND eq "Neg" ? "-" : $STRAND eq "Unk" ? "+" : $STRAND;
	print "\n>Doing ${YW}gene$N=$gene $LPR$gtfFolder$N\n - ${LCY}callFile$N=$callFile\n - ${LCY}output$N=$output\n";
	open (my $callFileIn, "<", $callFile) or DIELOG($outLog, "Failed to read from $callFile: $!\n");
	my $linecount = 0;
	my %peak; my ($minborder0, $maxborder0) = (-1,-1);
	while (my $line = <$callFileIn>) {
		chomp($line); $linecount ++;
		my ($read, @vals) = split("\t", $line);
		my ($beg, $end) = (0,0);
		my $vals = join("", @vals);
#		001100
#		012345
#		2 to 4
		my ($id) = parse_readName($read, $outLog);
		#my ($num1, $num2, $num3) = $read =~ /^.*\.?m(\d+_\d+).+\/(\d+)\/(ccs|\d+_\d+)/;
		#DIELOG($outLog, "\n\nERROR AT PARSING NUMBERS from read=$LGN$read$N\n\n") if not defined $num1 or not defined $num2 or not defined $num3;
		#$num3 = "0" if $num3 eq "ccs";
		#my $num = "$num1$num2$num3";
      #$num =~ s/_//g;
		LOG($outLog, "\tRead=$read, id=$id\n") if $linecount <= 3;
#		print "\n!!$BU$read$N\n" if $id eq "1803291104111173190";

#		$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[9]$/;

		my ($seqborder0) = $vals =~ /^(0+)[1-9]/; $seqborder0 = defined $seqborder0 ? length($seqborder0) : 0;
		my ($seqborder1) = $vals =~ /^(.+[1-9])0+$/; $seqborder1 = defined $seqborder1 ? length($seqborder1) : length($vals);
		($beg, $end) = ($seqborder0 + $BEG, $seqborder0 + $BEG + 1);
		my ($beg0, $end0) = (-1,-1);#$seqborder0,$seqborder1);
		$minborder0 = $beg if $minborder0 eq -1 or $minborder0 > $beg;
		my $print = "$CHR\tNA\texon\t$beg\t$end\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";# if (keys %clust) == 0;
		for (my $i = $seqborder0; $i < $seqborder1; $i++) {
			my $val = $vals[$i];
			($beg, $end) = ($BEG+$i, $BEG+$i+1);
			$beg0 = $beg if $val =~ /[89]/ and $beg0 eq -1;
			$end0 = $end if $val =~ /[89]/;
			$print .= "$CHR\tNA\texon\t$beg\t$end\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n" if $val =~ /[89]/;
#			$print .= "$CHR\tNA\tCDS\t$beg\t$end\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n" if $val =~ /[67]/ if (keys %clust) == 0;
		#	print $out "$CHR\tNA\t5UTR\t$beg\t$end\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n" if $val =~ /[45]/;
		}
		($beg, $end) = ($seqborder1 + $BEG, $seqborder1 + $BEG + 1);
		$maxborder0 = $end if $maxborder0 eq -1 or $maxborder0 < $end;
#		print "BEG=$BEG, END=$END, maxborder = $maxborder0, end=$end\n" if $linecount % 100 == 0;
		my $mid0 = int(($beg0+$end0)/2);
		$print .= "$CHR\tNA\texon\t$beg\t$end\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";# if (keys %clust) == 0;
		if ((keys %clust) == 0) {
			$peak{$mid0}{$beg}{$end}{$read}{print} = $print;
			$peak{$mid0}{$beg}{$end}{$read}{id} = $id;
		}
		else {
			my $y = $clust{$id}; die "Undef y atid=$id\n" if not defined $y;
			push(@{$peak{$y}{print}}, $print);
			push(@{$peak{$y}{id}}, $id);
#			push(@{$peak{$y}{end}}, $seqborder1 + $BEG + $y);
		}
	}
	LOG($outLog, "\n\n" . date() . " Printing to output file $LCY$output$N\n\n");
	open (my $out, ">", $output) or DIELOG($outLog, "Failed to write to $output: $!\n");
	print $out "$header\n";
	my $pos = $END;
	if ((keys %clust) != 0) {
		foreach my $y (sort {$a <=> $b} keys %peak) {
			my $printarr = $peak{$y}{print}; die "died at y=$y not defined arr\n" if not defined $printarr;
			my $idarr = $peak{$y}{id}; die "died at y=$y not defined arr\n" if not defined $idarr;
#			my $endarr = $peak{$y}{end}; die "died at y=$y not defined arr\n" if not defined $endarr;
			for (my $i = 0; $i < @{$printarr}; $i++) {
				my $print = $printarr->[$i];
				my $id = $idarr->[$i];
				my $end0 = $maxborder0 + $y; #+ $BEG;
				my $end1 = $end0 + 1;
#				print "end0=$end0 = maxborder0=$maxborder0 + BEG=$BEG = y=$y\n" if $i % 100 == 0;#if $id eq "113867" and $gene eq "RPS24" and $STRAND eq "+" and $label eq "PCB7";#output =~ PCB7_geneRPS24_Pos_20_0.35_CH.PEAK";
				print $out "$CHR\tNA\texon\t$minborder0\t$minborder0\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";
				print $out "$print";
				print $out "$CHR\tNA\texon\t$end0\t$end1\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";
			}
		}
	}
	elsif ((keys %clust) == 0) {
		foreach my $mid (sort {$a <=> $b} keys %peak) {
			foreach my $beg (sort {$a <=> $b} keys %{$peak{$mid}}) {
				foreach my $end (sort {$a <=> $b} keys %{$peak{$mid}{$beg}}) {
					foreach my $read (sort keys %{$peak{$mid}{$beg}{$end}}) {
						my $print = $peak{$mid}{$beg}{$end}{$read}{print};
						my $id = $peak{$mid}{$beg}{$end}{$read}{id};
						my $currpos = $linecount < 500 ? $pos : int($END + ($pos-$END)/$linecount*500);
						#die "currpos = liencount=$linecount < 500 ? pos=$pos or END=$END + (pos=$pos - END=$END) / lnecount=$linecount * 500)\n" if $id eq "113867" and $outName eq "PCB7_geneRPS24_Pos_20_0.35_CH.PEAK";
						print $out "$CHR\tNA\texon\t$minborder0\t$minborder0\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";
						print $out "$print";
						print $out "$CHR\tNA\texon\t$currpos\t$currpos\t0\t$STRAND\t0\tgene_id \"$id\"; transcript_id \"$id\"\n";
						$pos ++;
					}
				}
			}
		}
	}
	close $out;
	close $callFileIn;
	system("gzip -f $output");
}

sub parse_GTF_footPeak_logFile {
	my ($footPeakFolder, $footPeak_logFile, $outLog) = @_;
	my %opts; my $debugprint = "\n\n$YW<--------- 0. PARSING LOGFILES ----------$N\n\n"; my $optprint = "";
	open (my $footPeak_logFileIn, "<", $footPeak_logFile) or DIELOG($outLog, "Cannot read from $footPeak_logFile: $!\n");
	while (my $line = <$footPeak_logFileIn>) {
		chomp($line);
		my $check = 0;
		if ($line =~ />Run Params$/) {
			$check ++;
			die "0: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 1;
			while ($line !~ />Options:$/) {
				my ($param, $value) = $line =~ /^([\w ]+[a-zA-Z]+)[ \t]+:[ \t]+(.+)$/;
				$line =~ />Run Params$/ ? $debugprint .= "\n$LCY 0.$N footPeak $line:\n" : (defined $param and $param !~ /^Run script/) ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
				$opts{footPeak}{$param} = $value if (defined $param and defined $value and $param !~ /^Run script/);
				if (defined $param and $param eq "Run script short") {
					my @values = split(" -", $value); shift(@values);
					foreach my $values (@values) {
						my ($param2, $value2) = $values =~ /^(\w) (.+)$/;
						$debugprint .= "footLoopFolder = $LCY$value2$N\n" if $values =~ /^n /;
						$opts{footLoop2}{n} = $value2 if $values =~ /^n /;
						$opts{footPeak}{$param2} = $value2;
					}
				}
				$line = <$footPeak_logFileIn>; chomp($line);
			}
			$optprint .= ">footPeak\n";
			foreach my $param (sort keys %{$opts{footPeak}}) {
				$optprint .= "$param=$opts{footPeak}{$param}\n";
			}
		}
		if ($line =~ /^>Options/) {
			$check ++;
			die "1: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 2;
			while ($line !~ />Run Params from footLoop.pl/) {
				my ($param, $desc, $value) = $line =~ /^\-(\w)\s*(\w+)\s*:\s*([a-zA-Z0-9]+)$/;
				   ($param, $value) = $line =~ /^\-(\w)\s*:\s*([a-zA-Z0-9]+)$/ if not defined $param;
				$line =~ />Options/ ? $debugprint .= "\n$LGN 1.$N footPeak $line: " : defined $param ? $debugprint .= "$LCY$param$N=$value;" : $debugprint .= "";
				$opts{footPeak2}{$param} = $value if (defined $param and defined $value);
				$line = <$footPeak_logFileIn>; chomp($line);
			}
			foreach my $param (sort keys %{$opts{footPeak2}}) {
				$optprint .= "$param=$opts{footPeak2}{$param}\n";
			}
		}
		if ($line =~ />Run Params from footLoop.pl/) {
			$check ++; my $debugprint2;
			die "2: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 3;
			while ($line !~ /^>Options from footLoop.pl/) {
				my ($param, $value) = $line =~ /^footLoop\s*(\w+)\s*:[ \t]+(.+)$/;
				if ($line =~ /from footLoop.pl logfile/) {($param, $value) = $line =~ /(logfile)=\.?\/?(.+)$/; $value = $opts{footLoop2}{n} . $value; $debugprint2 .= "\t- $LCY$param$N=$value\n";}
				$line =~ />Run Params/ ? $debugprint .= "\n$LCY 2.$N footLoop $line:\n$debugprint2" : defined $param ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
				$opts{footLoop}{$param} = $value if (defined $param and defined $value and $param ne "origDir");
				$line = <$footPeak_logFileIn>; chomp($line);
			}
			$optprint .= ">footLoop\n";
			foreach my $param (sort keys %{$opts{footLoop}}) {
				$optprint .= "$param=$opts{footLoop}{$param}\n";
			}

		}
		if ($line =~ /^>Options from footLoop.pl/) {
			$check ++;
			die "3: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 4;
			while ($line !~ /^.+[\-]+\>/) {
				my ($param, $value) = $line =~ /^\-(\w)\s*:\s+(.+)$/;
				$opts{footLoop2}{$param} = $value if (defined $param and defined $value);
				$line =~ />Options from footLoop/ ? $debugprint .= "\n$YW 3.$N footLoop $line:\n" : defined $param ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
				$line = <$footPeak_logFileIn>; chomp($line);
			}
			foreach my $param (sort keys %{$opts{footLoop2}}) {
				$optprint .= "$param=$opts{footLoop2}{$param}\n";
			}
		}
	#	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $bedLine);
	}
	$debugprint .= "\n\n$YW------------------------------------->$N\n\n\n";
	makedir("$footPeakFolder/GTF");
	open (my $out, ">", "$footPeakFolder/GTF/.PARAMS") or DIELOG($outLog, "Failed to write to $footPeakFolder/GTF/.PARAMS: $!\n");
	print $out "$optprint\n";
	close $out;
	LOG($outLog, $debugprint);
	return \%opts;
}

sub parse_geneIndexFile {
	my ($footLoopFolder, $outLog) = @_;
	my ($geneIndexFile) = <$footLoopFolder/*.bed>;
	my %coor;
	LOG($outLog, "${LCY}geneIndexFile$N=$geneIndexFile\n");
	die "geneindexFile does not exist!\n" if not defined $geneIndexFile;
	open (my $in, "<", $geneIndexFile) or DIELOG($outLog, "Failed to read from $geneIndexFile: $!\n");
	while (my $line = <$in>) {
		chomp($line);
		my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
		$gene = uc($gene);
		$coor{$gene}{chr} = $chr;
		$coor{$gene}{beg} = $beg;
		$coor{$gene}{end} = $end;
		$coor{$gene}{strand} = $strand;
		#print "chr=$chr:$beg-$end gene=$gene, strand=$strand\n";
	}
	close $in;
	return \%coor;
}

sub parse_fasta {
	my %data;
	my ($output) = @_;
	# get fasta
	my $fasta = $output . "fa";my $fastaCount = 0;
	while (-e $fasta) {
		$fastaCount ++;
		$fasta = $output . $fastaCount . "fa";
	}
	#not stranded because we agree to use relative to + strand
	my %head;
	my $facmd = "fastaFromBed -fi $opt_f -bed $opt_i -fo $fasta -name";
	system($facmd) == 0 or die "Failed to run $CY$facmd$N: $!\n";
	open (my $faIn, "<", $fasta) or die "Failed to read from $fasta: $!\n";
	$fasta = new FAlite($faIn);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;$def =~ s/^>//;my $gene = $def;
		my @seq = split("", uc($entry->seq));
		my $beg = $data{$gene}{beg};die if not defined $beg;
		my ($beg0, $beg1) = ($beg, $beg + 1);#($beg - 500, $beg - 500+ 1);$beg0 = 1 if $beg0 < 0;$beg1 = 1 if $beg1 < 1;
		my $end = $data{$gene}{end};die if not defined $end;
		my ($end0, $end1) = ($end - 1, $end); #($end - 500, $end - 500+ 1);$end0 = 1 if $end0 < 0;$end1 = 1 if $end1 < 1;
		my $chr = $data{$gene}{chr};die if not defined $chr;
		my $strand = $data{$gene}{strand};die if not defined $strand;
		$head{BLOCK}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.BLOCK\"; transcript_id \"MYNUMBER.$gene.BLOCK\"\n";
		$head{BLOCK}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.BLOCK\"; transcript_id \"MYNUMBER.$gene.BLOCK\"\n";
		$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		for (my $i = 0;$i < @seq;$i++) {
			my $begPeak = $i + $beg;# - $bufferL;# 0 based
			my $endPeak = $begPeak + 1;# 1 based
			if ($seq[$i] eq "C" and defined $seq[$i+1] and $seq[$i+1] eq "G") {
				$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t" . ($begPeak+1) . "\t" . ($endPeak+1) . "\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
			}
			elsif ($seq[$i] eq "C") {
				$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
			}
			elsif ($seq[$i] eq "G" and $i != 0 and defined $seq[$i-1] and $seq[$i-1] ne "C") {
				$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
			}
		}
		$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
	}
#	foreach my $type (sort keys %head) {
#		next if $type eq "BLOCK";
#		foreach my $strand (sort keys %{$head{$type}}) {
#			print $GTFPosfile "track name=\"$type\_header +\" itemRgb=On\n" if $strand eq "pos";
#			print $GTFPosfile "track name=\"$type\_header +\" itemRgb=On\n" if $strand eq "pos";
#			print $GTFNegfile "track name=\"$type\_header -\" itemRgb=On\n" if $strand eq "neg";
#			print $GTFNegfile "track name=\"$type\_header -\" itemRgb=On\n" if $strand eq "neg";
#			foreach my $gene (sort keys %{$head{$type}{$strand}}) {
#				print $GTFPosfile "$head{$type}{$strand}{$gene}\n" if $strand eq "pos";
#				print $GTFNegfile "$head{$type}{$strand}{$gene}\n" if $strand eq "neg";
#			}
#		}
#	}
	close $faIn;
	return(\%head);
}


# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

sub parse {
	my %head;
	my ($gene, $file, $coor, $out) = @_;
	print "Parsing gene=$gene, file=$file\n";
	my $chr = $coor->{chr};
	my $beg = $coor->{beg};my ($beg0, $beg1) = ($beg, $beg + 1);
	my $end = $coor->{end};my ($end0, $end1) = ($end - 1, $end);
	my $strand = $coor->{strand};
	my $peak = 0;
	#hclust_R($file);# produces $file.Clust
	my %peak;
	open (my $in, "<", "$file") or die "Failed to read from $file.Clust: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		$line =~ s/["']//g;$line =~ s/^SEQ_//;
		my ($name, @vals) = split ("\t", $line);
		$name = "MYNUMBER.$gene.$name";
		#my $begGTF;# = "$chr\tNA\texon\t$beg0\t$beg1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		#my $endGTF;# = "$chr\tNA\texon\t$end0\t$end1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		my $begGTF = "$chr\tNA\texon\t$beg0\t$beg1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		my $endGTF = "$chr\tNA\texon\t$end0\t$end1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		$peak{$peak}{print} .= "$begGTF\n";
		my ($peakMinI, $peakMaxI, $peakMin, $peakMax) = (int($end/200+0.5), int($beg/200+0.5), $end, $beg);#to determine cluster
		for (my $i = 0;$i < @vals;$i++) {
			next if $vals[$i] =~ /^[62]$/;
			my $begPeak = $i + $beg + $bufferL;# 0 based
			my $endPeak = $begPeak + 1;# 1 based

			# first and last coordinate where the read has data
			$begGTF = $begPeak if not defined $begGTF;
			$endGTF = $endPeak;

			$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[9]$/;
#			$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[1]$/;
			$peak{$peak}{print} .= "$chr\tNA\tCDS\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[5]$/;
			#$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[4]$/;
			#$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[0]$/;
			
			# For sorting
			my $begPeakIndex = int($begPeak/200+0.5);
			my $endPeakIndex = int($endPeak/200+0.5);
			if ($peakMinI > $begPeakIndex and $vals[$i] =~ /^[59]$/) {
				$peakMinI = $begPeakIndex;
				$peakMin  = $begPeak;
			}
			if ($peakMaxI < $endPeakIndex and $vals[$i] =~ /^[59]$/) {
				$peakMaxI = $endPeakIndex;
				$peakMax  = $endPeak;
			}
		}
		#my $begGTF0 = $begGTF; my $begGTF1 = $begGTF + 1;
		#my $endGTF0 = $endGTF - 1; my $endGTF1 = $endGTF;
		#$peak{$peak}{print} = "$chr\tNA\texon\t$begGTF0\t$begGTF1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" . $peak{$peak}{print};
		$peak{$peak}{print} .= "$endGTF\n";
		#$peak{$peak}{print} .= "$chr\tNA\texon\t$endGTF0\t$endGTF1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n";
		$peak{$peak}{begI} = $peakMinI;
		$peak{$peak}{endI} = $peakMaxI;
		$peak{$peak}{beg} = $peakMin;
		$peak{$peak}{end} = $peakMax;
		$peak++;
	}
	
	my $countz = 0;
	my $strand2 = $strand eq "+" ? "pos" : $strand eq "-" ? "neg" : die "Cannot determine strand ($strand), is not + or - !\n";
	foreach my $peak (sort {$peak{$b}{begI} <=> $peak{$a}{begI} || $peak{$b}{endI} <=> $peak{$a}{endI} || $peak{$b}{beg} <=> $peak{$a}{beg} || $peak{$b}{end} <=> $peak{$a}{end}} keys %peak) {
		my $countz2 = int($countz/100);
		my $PEAK = $peak{$peak}{print}; $PEAK =~ s/MYNUMBER/$countz2/g;
		print $out $PEAK;
		if ($countz % 100 == 0) {
			die "type = CHH, GENE = $gene, STRAND=$strand2\n" if not defined $head{CHH}{$strand2}{$gene};
			my $CHH = $head{CHH}{$strand2}{$gene}; $CHH =~ s/MYNUMBER/$countz2/g;
			my $CPG = $head{CPG}{$strand2}{$gene}; $CPG =~ s/MYNUMBER/$countz2/g;
			my $BLOCK = $head{BLOCK}{$strand2}{$gene}; $BLOCK =~ s/MYNUMBER/$countz2/g;
			print $out $CHH;
			print $out $CPG;
			print $out $BLOCK;
		}
		$countz ++;
	}
}

sub hclust_R {
	my ($reads) = @_;
	my $output = "$reads.Clust";
	open (my $out, ">", "$output.R");
	print $out "
	library(\"GMD\")
	df = read.table(\"./$reads\", sep=\"\t\", row.names=1)
	df2 = df
	df2[df2 != 9] = 0
	df\$sum = apply(df2,1,sum)
	df = df[order(-df\$sum),]
	dimz = dim(df)[1]
	df = df[1:dimz,]
	df = subset(df,select=-sum)
	df2 = df
	df2[df2 !=9] = 0
	h = heatmap.3(
		df2,
		dendrogram=\"row\",
		Rowv=TRUE, Colv=FALSE,
		labRow=FALSE, labCol=FALSE
	)
	id = h\$rowInd
	id = rev(id)
	df = df[id,]
	dev.off()
	write.table(df, file = \"./$output\", sep = \"\t\")
	";
	close $out;
	system("R --vanilla --no-save < $output.R") == 0 or die "Failed to run Rscript hclust ($output.R): $!\n";
	system("rm $output.R") == 0 or die "Failed to remove $output.R: $!\n";
}
