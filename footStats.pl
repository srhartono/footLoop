#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G);
getopts("vd:n:G:");

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
#  print "\n\n\e[1;33m ------------ BEGIN ------------> \e[0m\n";
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
my $date = getDate();
my $uuid = getuuid();
my $numThread = 1;
my ($footPeakFolder, $outLog) = check_sanity();
my $genewant = $opt_G;
my $opts = parse_footPeak_logFile($footPeakFolder, "$footPeakFolder/footPeak_logFile.txt", $outLog);
my ($footLoopFolder) = $opts->{footLoop2}{n};
my $geneIndexFile = $opts->{footLoop2}{i};
my $coor = parse_geneIndexFile($geneIndexFile, $outLog);
DIELOG($outLog, date() . " ERROR: $footLoopFolder footloop folder doesn't exist!\n") if not -d $footLoopFolder;

# get sam file
my ($samFile) = $opts->{footLoop}{samFile};
my (%samData) = %{parse_samFile($samFile, $outLog)};
foreach my $gene (sort keys %samData) {
	my @key = qw(Pos Neg);
	my $total = $samData{$gene}{total};
	my $pos = $samData{$gene}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $samData{$gene}{Neg}; $neg = 0 if not defined $neg;
#	print "$gene\t$total\t$pos\t$neg\n";
}
my ($samFixedFile) = $opts->{footLoop}{samFixed};
my ($origDir) = getFilename($samFixedFile, "folder");
my (%samFixedData) = %{parse_samFile($samFixedFile, $outLog)};

my %print0;
$print0{header} = "gene\tsam_uread\tsam_uread_pos\tsam_uread_neg\tsamfix_uread_pos\tsamfix_uread_neg\tsamfix_pos_pos\tsamfix_neg_neg\tsamfix_pos_neg\tsamfix_neg_pos\tsamfinal_uread_pos\tsamfinal_uread_neg\tsamfinal_uread_unk";
foreach my $gene (sort keys %samFixedData) {
	my @key = qw(Pos Neg);
	my $total = $samFixedData{$gene}{total};
	my $sampos = $samData{$gene}{Pos}; $sampos = 0 if not defined $sampos;
	my $samneg = $samData{$gene}{Neg}; $samneg = 0 if not defined $samneg;
	my $pos = $samFixedData{$gene}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $samFixedData{$gene}{Neg}; $neg = 0 if not defined $neg;
	die "Discrepancy betweensam and fixed sam pos and neg in gene=$gene (pos=$pos, sampos=$sampos)\n" if $sampos ne $pos;
	die "Discrepancy betweensam and fixed sam neg and neg in gene=$gene (neg=$neg, samneg=$samneg)\n" if $samneg ne $neg;
	my $pos2 = $samFixedData{$gene}{strand}{Pos}; $pos2 = 0 if not defined $pos2;
	my $neg2 = $samFixedData{$gene}{strand}{Neg}; $neg2 = 0 if not defined $neg2;
	my $posneg = $samFixedData{$gene}{change}{Pos}{Neg}; $posneg = 0 if not defined $posneg;
	my $negpos = $samFixedData{$gene}{change}{Neg}{Pos}; $negpos = 0 if not defined $negpos;
	my $pospos = $samFixedData{$gene}{change}{Pos}{Pos}; $pospos = 0 if not defined $pospos;
	my $negneg = $samFixedData{$gene}{change}{Neg}{Neg}; $negneg = 0 if not defined $negneg;
	my $posorig = "$origDir/$gene\_Pos.orig"; ($posorig) = -e $posorig ? `wc -l $posorig` =~ /^(\d+)/ : 0;
	my $negorig = "$origDir/$gene\_Neg.orig"; ($negorig) = -e $negorig ? `wc -l $negorig` =~ /^(\d+)/ : 0;
	my $unkorig = "$origDir/$gene\_Unk.orig"; ($unkorig) = -e $unkorig ? `wc -l $unkorig` =~ /^(\d+)/ : 0;
	$print0{gene}{$gene} = "$gene\t$total\t$pos\t$neg\t$pos2\t$neg2\t$pospos\t$negneg\t$posneg\t$negpos\t$posorig\t$negorig\t$unkorig";
}
sub parse_samFile {
	my ($samFile, $outLog) = @_;
	my %data;
	my $in;
	if ($samFile =~ /.gz$/) {
		open ($in, "zcat $samFile|") or DIELOG($outLog, date() . " Failed to read from samFile $LCY$samFile$N: $!\n");
	}
	else {
		open ($in, "<", "$samFile") or DIELOG($outLog, date() . " Failed to read from samFile $LCY$samFile$N: $!\n");
	}
	while (my $line = <$in>) {
		chomp($line);
		my @arr = split("\t", $line);
		next if @arr < 6;
		my ($read, $flag, $gene, $pos, $qual) = @arr;
		my ($flag2, $type);
		if ($arr[1] !~ /^(0|16)$/) {
			($read, $type, $flag, $flag2, $gene) = @arr;
		}
		my $strand = $flag eq 16 ? "Neg" : $flag eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$samFile$N, flag=$flag, line:\n\n$line\n\n");
		$data{$gene}{$strand} ++;
		$data{$gene}{total} ++;
		if ($arr[1] !~ /^(0|16)$/) {
			my $strand2 = $flag2 eq 16 ? "Neg" : $flag2 eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$samFile$N, flag2=$flag2, line:\n\n$line\n\n");
			$data{$gene}{strand}{$strand2} ++;
			$data{$gene}{change}{$strand}{$strand2} ++;
		}
	}
	close $in;
	return \%data;
}

# get fixed sam file
my ($fixedsamFile) = $opts->{footLoop}{samFixed};



my (@TXTFile) = <$footPeakFolder/99_FOOTSTATS/.PEAKSTATSTEMP/.0_*.TXT>;
if (@TXTFile == 0) {@TXTFile = <$footPeakFolder/.0_*.TXT>;}
open (my $out, ">", "$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT") or DIELOG($outLog, date() . " Failed to write to $footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT: $!\n");

my $data;
foreach my $TXTFile (@TXTFile) {
	my ($folder1, $fileName1) = getFilename($TXTFile, "folderfull");
	my ($footName) = $fileName1 =~ /.0_RESULTS_(.+).TXT$/;
	print "Parsed file=$fileName1 footName = $LCY$footName$N\n";
	my @result = parseName($footName . "_CH");
	my @header = qw(label gene strand window thres conv bc plasmid desc);
	my ($label, $gene, $strand, $window, $thres, $conv, $bc, $plasmid, $desc) = @result;
	my $label2 = "$label\_$bc\_$plasmid\_$desc" if defined $bc and $bc ne ""; 
	$label2 = $label if not defined $bc or (defined $bc and $bc eq "");
	$result[5] = "";
	$data = parseTXT($data, $TXTFile, $label2, $gene, $strand, $window, $thres);
}

my %print1;
my %print;
my $LABEL;
print $out "type\tlabel\tgene\tread_strand\twindow\tthres\tconv\tread_unique_total\tread_unique_peak_total\tread_unique_peak_perc\tfootPeakFolder\tpeakFile\n";
foreach my $label (sort keys %{$data}) {
	foreach my $gene (sort keys %{$data->{$label}}) {
		foreach my $read_strand (sort keys %{$data->{$label}{$gene}}) {
			foreach my $window (sort keys %{$data->{$label}{$gene}{$read_strand}}) {
				foreach my $thres (sort keys %{$data->{$label}{$gene}{$read_strand}{$window}}) {
					foreach my $conv (sort keys %{$data->{$label}{$gene}{$read_strand}{$window}{$thres}}) {
						my $strand = $coor->{uc($gene)}{strand}; $strand = "Pos" if $strand eq "+"; $strand = "Neg" if $strand eq "-";
						my $goodconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "CG" : $read_strand eq "Neg" ? "GC" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						$goodconv = $goodconv ne "H" ? $goodconv : $read_strand eq "Pos" ? "CH" : $read_strand eq "Neg" ? "GH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						my $revconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "GC" : $read_strand eq "Neg" ? "CG" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						$revconv = $revconv ne "H" ? $revconv : $read_strand eq "Pos" ? "GH" : $read_strand eq "Neg" ? "CH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						my $type;
						$type = 	$read_strand eq "UNK" ? "UNK" :
									($strand eq $read_strand and $conv eq $goodconv) ? "PEAK" :
									($strand ne $read_strand and $conv eq $goodconv) ? "PEAKNEG" : 
									($strand eq $read_strand and $conv eq $revconv) ? "PEAKREV" :
									($strand ne $read_strand and $conv eq $revconv) ? "PEAKNEGREV" : "OTHERS";
						my ($currfootPeakFolder) = getFullpath($data->{$label}{$gene}{$read_strand}{$window}{$thres}{$conv}{footPeakFolder});
						my ($currpeakFile) = $data->{$label}{$gene}{$read_strand}{$window}{$thres}{$conv}{peakFile};
						my ($read_unique_total) = $data->{$label}{$gene}{$read_strand}{$window}{$thres}{$conv}{read_unique_total};
						my ($read_unique_peak_total) = $data->{$label}{$gene}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_total};
						my ($read_unique_peak_perc) = $data->{$label}{$gene}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_perc};
						$print1{$gene}{$type} += $read_unique_peak_total;
						$print1{$gene}{strand} = $strand;
						$print1{$gene}{OTHERS_TOTAL} ++ if $type eq "OTHERS";
						$print{$type} .= "$type\t$label\t$gene\t$read_strand\t$window\t$thres\t$conv\t$read_unique_total\t$read_unique_peak_total\t$read_unique_peak_perc\t$currfootPeakFolder\t$currpeakFile\n";
						$LABEL = $label;
					}
				}
			}
		}
	}
}

open (my $out2, ">", "$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT") or die;
print $out2 "label\tstrand\t$print0{header}";
my @types = qw(PEAK PEAKNEG PEAKREV PEAKNEGREV OTHERS); 
for (my $i = 0; $i < @types; $i++) { 
	print $out $print{$types[$i]};
	print $out2 "\t$types[$i]";
}
print $out2 "\n";

foreach my $gene (sort keys %{$print0{gene}}) {
	my $print0 = $print0{gene}{$gene};
	my $strand = $coor->{uc($gene)}{strand}; $strand = "Pos" if $strand eq "+"; $strand = "Neg" if $strand eq "-";
	print $out2 "$LABEL\t$strand\t$print0";
	for (my $i = 0; $i < @types-1; $i++) { 
		my $print1 = $print1{$gene}{$types[$i]}; $print1 = 0 if not defined $print1;
		print $out2 "\t$print1";
	}
	my $printOTHERS = $print1{$gene}{OTHERS}; $printOTHERS = 0 if not defined $printOTHERS;
	my $printOTHERStot = $print1{$gene}{OTHERS_TOTAL}; $printOTHERStot = 1 if not defined $printOTHERStot;
	$printOTHERS = int($printOTHERS / $printOTHERStot * 10)/10;
	print $out2 "\t$printOTHERS\n";
}
LOG($outLog, date() . "$LCY$opt_n$N Done!\n\nOutput:
$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT
$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT\n");

sub check_sanity {
	die "\nUsage: $YW$0$N [Optional: -G <genewant>] -n $LCY$opt_n$N\n" unless defined $opt_n and -d $opt_n;
	my ($footPeakFolder) = $opt_n;
	my ($logFile) = $footPeakFolder . "/footStats_logFile.txt";
	my ($outDir) = $footPeakFolder . "/99_FOOTSTATS/";
	my $optG = defined $opt_G ? "-G $opt_G " : "";
	makedir("$footPeakFolder/99_FOOTSTATS/") if not -d "$footPeakFolder/99_FOOTSTATS/";
	open (my $outLog, ">", $logFile) or die "\nFailed to write to $LCY$logFile$N: $!\n";
	print $outLog "$0 version $version\n\nRun script: $0 $optG-n $opt_n\n";
	return($footPeakFolder, $outLog);
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

sub parseTXT {
	my ($data, $TXTFile, $label, $gene, $strand, $window, $thres) = @_;
	my @header = qw(footPeakFolder peakFile gene conv read_unique_total read_unique_peak_total read_unique_peak_perc);

	open (my $in1, "<", $TXTFile) or die "Cannot read from $TXTFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		if ($line =~ /^#/) {
			next;
		}
		my @arr = split("\t", $line);
		die "total number of header isn't the same as total column in file $LCY$TXTFile$N!\n\n$LPR$line$N\n\n" if @header != @arr;
#		print "\n$line\n";
		my $conv = $arr[3];
		for (my $i = 0; $i < @arr; $i++) {
			next if $i == 3 or $i == 2;
			$data->{$label}{$gene}{$strand}{$window}{$thres}{$conv}{$header[$i]} = $arr[$i];
		}
	}
	close $in1;
	return($data);
}

sub parse_footPeak_logFile {
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
        #       my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $bedLine);
        }
        $debugprint .= "\n\n$YW------------------------------------->$N\n\n\n";
        makedir("$footPeakFolder/99_FOOTSTATS");
        open (my $out, ">", "$footPeakFolder/99_FOOTSTATS/.PARAMS") or DIELOG($outLog, "Failed to write to $footPeakFolder/99_FOOTSTATS/.PARAMS: $!\n");
        print $out "$optprint\n";
        close $out;
        LOG($outLog, $debugprint);
        return \%opts;
}

sub parse_geneIndexFile {
   my ($geneIndexFile, $outLog) = @_;
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

__END__
my %data;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

170802_pcb05_RnaseHneg_ccs_3minFP.fastq.gz.rmdup.fq.gz_PEAK_20W35T      CALM3_Neg_20_0.35_GC.PEAK       CALM3   GC      160     3       1.9
