#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_c $opt_F);
getopts("vd:n:G:cF");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
}

use myFootLib;
use FAlite;

my $homedir = $ENV{"HOME"};

my ($footLoop_script_folder, $version, $md5script) = check_software();
my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
$version_small = $version if not defined $version_small;

#my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
#my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
#my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
#my $version = join("", @version);
#if (defined $opt_v) {
#   print "$version\n";
#   exit;
#}
#my ($version_small) = "vUNKNOWN";
#foreach my $versionz (@version[0..@version-1]) {
#   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
#}

my $date = getDate();
my $uuid = getuuid();
my $numThread = 1;
my ($footPeakFolder, $outLog) = check_sanity();
my $genewant = $opt_G;
my $opts = footStats_parse_footPeak_logFile($footPeakFolder, "$footPeakFolder/footPeak_logFile.txt", $outLog);
my ($footLoopFolder) = $opts->{footLoop2}{n};
my $geneIndexFile = $opts->{footLoop2}{i};
my $coor = parse_geneIndexFile($geneIndexFile, $outLog);
DIELOG($outLog, date() . " ERROR: $footLoopFolder footloop folder doesn't exist!\n") if not -d $footLoopFolder;
my %data;
my ($labelFile) = "$footPeakFolder/.LABEL";

#for (my $i = 0; $i < @peakFiles; $i++) {
#	my $peakFile = $peakFiles[$i];
#	print "$i $peakFile\n";
#	my $file = $peakFile;
#	print "file=$file\n";
#	open (my $in1, "<", $file) or LOG($outLog, "Failed to read from $file: $!\n") and exit 1;
#	while (my $line = <$in1>) {
#		chomp($line);
#		my ($read, @val) = split("\t", $line);
#		print "HERE\n$read\n";
#		last;
#	}
#	close $in1;
#	last if $i > 10;
#}
my $wclFile = "$footPeakFolder/99_FOOTSTATS/0_wcl.tsv";
if (not -e $wclFile or defined $opt_F) {
	my $cmd = "wc -l $footPeakFolder/.CALL/*.out > $wclFile";
	print "\nGetting linecount of each .out files:\n$LCY$cmd$N\n\n";
	system($cmd) == 0 or die "Failed to run cmd: $!\n$LCY$cmd$N\n\n";
}
my $wclhash = parse_wcl($wclFile);
my %wcl = %{$wclhash};
sub parse_wcl {
	my ($wclFile) = @_;
	my %wcl;
	open (my $in1, "<", $wclFile) or die "Failed to read from $wclFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		$line =~ s/^\s+(\d+)/$1/;
		my ($totalread, $file) = split(" ", $line);
		#print "$file\t$totalread\n";
		$wcl{$file} = $totalread;
	}
	close $in1;
	return(\%wcl);
}

my @peakFiles = <$footPeakFolder/.CALL/*PEAK*out>;
my %final;
my ($len) = $peakFiles[0] =~ /^.+_l(\d+)/;
for (my $i = 0; $i < @peakFiles; $i++) {
	print "Done $i\n" if $i % 1000 == 0;
	my $peakFile = $peakFiles[$i];
	($data{$peakFile}{total}) = defined($wcl{$peakFile}) ? $wcl{$peakFile} : 0;
	#($data{$peakFile}{total}) = `wc -l $peakFile` =~ /^(\d+)/ if -e $peakFile;
	my $parseName = parseName($peakFile);
	my ($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	my $readstrand = $strand;
	$label =~ s/_bismark_bt2.bam//;
	my $plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	my $genestrand = $coor->{$plasmid2}{strand}; 
	$genestrand = "Pos" if not defined $genestrand;
	$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	my $array = $data{$peakFile}{array};
   my ($flag) = getFlag($peakFile, $genestrand, $readstrand, $convtype);
	#print "$i $YW$bc $plasmid2 $plasmid $desc $LPR$genestrand $LCY$readstrand $convtype $LPR$flag$N\t$LGN$data{$peakFile}{total}$N\t$peakFile\n";

	my $nopkFile = $peakFile; $nopkFile =~ s/PEAK/NOPK/g;
	($data{$nopkFile}{total}) = defined($wcl{$nopkFile}) ? $wcl{$nopkFile} : 0;
	#($data{$nopkFile}{total}) = `wc -l $nopkFile` =~ /^(\d+)/ if -e $nopkFile;
	my $totalread = $data{$peakFile}{total} + $data{$nopkFile}{total};
	my $print = "$bc\t$plasmid\t$desc\t$readstrand\t$convtype\t$data{$peakFile}{total}\t$totalread\t$flag\t$peakFile\n";
	my $plasmiddesc = "PLASMID$plasmid;DESC$desc";
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{flag} = $flag;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak} = $data{$peakFile}{total};
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk} = $data{$nopkFile}{total};
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakperc} = $totalread == 0 ? 0 : int(1000*$data{$peakFile}{total}/$totalread+0.5)/10;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkperc} = $totalread == 0 ? 0 : int(1000*$data{$nopkFile}{total}/$totalread+0.5)/10;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread} = $totalread;
	#$parseName = parseName($nopkFile);
	#($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	#$readstrand = $strand;
	#$label =~ s/_bismark_bt2.bam//;
	#$plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	#$genestrand = $coor->{$plasmid2}{strand}; 
	#$genestrand = "Pos" if not defined $genestrand;
	#$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	#$array = $data{$nopkFile}{array};
   #($flag) = getFlag($nopkFile, $genestrand, $readstrand, $convtype);
	#print "$i\t$YW$bc\t$plasmid\t$desc\t$readstrand\t$convtype\t$LGN$data{$peakFile}{total}$N\t$LGN$totalread$N\t$LPR$flag$N\n";
	#print "$i $YW$bc $plasmid2 $plasmid $desc $LPR$genestrand $LCY$readstrand $convtype $LPR$flag$N\t$LGN$data{$nopkFile}{total}$N\t$nopkFile\n\n";
}

my @strand = qw(Pos Neg);
my @convtype = qw(CG CH GC GH);
my @strandconvtype = qw(Pos_CG Neg_GC Pos_GC Neg_CG Pos_CH Neg_GH Pos_GH Neg_CH);
open (my $out1, ">", "$footPeakFolder/99_FOOTSTATS/1_STATS.txt") or die;
print $out1 "BC\tplasmid\tdesc\ttx";
#foreach my $strand (@strandtrand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
foreach my $strand (@strand[0..@strand-1]) {
	print $out1 "\t$strand";
}
#	}
#}
#foreach my $strand (@strand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\t$strandconvtype";
}
#	}
#}
print $out1 "\n";
foreach my $bc (sort keys %final) {
	foreach my $plasmiddesc (sort keys %{$final{$bc}}) {
		my ($plasmid, $desc) = $plasmiddesc =~ /^PLASMID(.+);DESC(.+)$/;
		my ($tx) = $desc =~ /(TX_RH1|TX|NT|NT_RH1)$/;
		my ($desc2) = $desc =~ /^(.+)_(TX_RH1|TX|NT|NT_RH1)$/;
		$tx = "NA" if not defined $tx;
		$desc2 = $desc if not defined $desc2;
		print $out1 "$bc\t$plasmid\t$desc2\t$tx";
		my $peakpercprint = "";
		my $totalreadprint = "";
		foreach my $strand (@strand[0..@strand-1]) {	
			my $convtype = "CH";
			my $totalread = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread};
			$totalread = 0 if not defined $totalread;
			$totalreadprint .= "\t$totalread";
		}
		foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
		#foreach my $strand (@strand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
		#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
				my ($strand, $convtype) = $strandconvtype =~ /^(\w+)_(\w+)$/;
				my $flag = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{flag};				
				my $peak = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak};				
				my $nopk = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk};				
				my $totalread = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread};				
				my $peakperc = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakperc};				
				my $nopkperc = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkperc};
				$totalread = 0 if not defined $totalread;
				$peakperc = 0 if not defined $peakperc;
				$nopkperc = 0 if not defined $nopkperc;
				$peakpercprint .= "\t$peakperc";
		#	}
		#}
		}

		print $out1 "$totalreadprint";
		print $out1 "$peakpercprint\n";
	}
}

print "SUCCESS!!!\n";
exit 0;
sub footStats_parse_footPeak_logFile {
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
	#	   my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $bedLine);
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
	  $coor{uc($gene)}{chr} = $chr;
	  $coor{uc($gene)}{beg} = $beg;
	  $coor{uc($gene)}{end} = $end;
	  $coor{uc($gene)}{strand} = $strand;
	  #print "chr=$chr:$beg-$end gene=$gene, strand=$strand\n";
   }
   close $in;
   return \%coor;
}

sub check_sanity {

	my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N [Optional: -G <genewant>] -n $LCY<footPeak_folder>$N

";

	die $usage unless defined $opt_n and -d $opt_n;

	my ($footPeakFolder) = $opt_n;
	my ($logFile) = $footPeakFolder . "/footStats_logFile.txt";
	my ($outDir) = $footPeakFolder . "/99_FOOTSTATS/";
	my $optG = defined $opt_G ? "-G $opt_G " : "";
	makedir("$footPeakFolder/99_FOOTSTATS/") if not -d "$footPeakFolder/99_FOOTSTATS/";
	open (my $outLog, ">", $logFile) or die "\nFailed to write to $LCY$logFile$N: $!\n";
	print $outLog "$0 version $version\n\nRun script: $0 $optG-n $opt_n\n";
	return($footPeakFolder, $outLog);
}

=comment
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Processing BAMFile$N $LCY$BAMFile$N\n");
my (%BAMData) = %{parse_BAMFile($BAMFile, $outLog)};
foreach my $gene (sort keys %BAMData) {
	my @key = qw(Pos Neg);
	my $total = $BAMData{uc($gene)}{total};
	my $pos = $BAMData{uc($gene)}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $BAMData{uc($gene)}{Neg}; $neg = 0 if not defined $neg;
#	print "$gene\t$total\t$pos\t$neg\n";
}

my ($BAMFixedFile) = $opts->{footLoop}{BAMFixed};
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Processing BamFixedFile$N $LCY$BAMFixedFile$N\n");
my ($origDir) = getFilename($BAMFixedFile, "folder");
my (%BAMFixedData) = %{parse_BAMFile($BAMFixedFile, $outLog)};

LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Printing each gene in BamFixedData$N\n");
my %dataorig;
my %print0;
$print0{header} = "gene\tsam_uread\tsam_uread_pos\tsam_uread_neg\tsamfix_uread_pos\tsamfix_uread_neg\tsamfix_pos_pos\tsamfix_neg_neg\tsamfix_pos_neg\tsamfix_neg_pos\tsamfinal_uread_pos\tsamfinal_uread_neg\tsamfinal_uread_unk\tinclude_cpg";
foreach my $gene (sort keys %BAMFixedData) {
	my @key = qw(Pos Neg);
	my $total = $BAMFixedData{uc($gene)}{total};
	my $sampos = $BAMData{uc($gene)}{Pos}; $sampos = 0 if not defined $sampos;
	my $samneg = $BAMData{uc($gene)}{Neg}; $samneg = 0 if not defined $samneg;
	my $pos = $BAMFixedData{uc($gene)}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $BAMFixedData{uc($gene)}{Neg}; $neg = 0 if not defined $neg;
	die "Discrepancy betweensam and fixed sam pos and neg in gene=$gene (pos=$pos, sampos=$sampos)\n" if $sampos ne $pos;
	die "Discrepancy betweensam and fixed sam neg and neg in gene=$gene (neg=$neg, samneg=$samneg)\n" if $samneg ne $neg;
	my $pos2 = $BAMFixedData{uc($gene)}{strand}{Pos}; $pos2 = 0 if not defined $pos2;
	my $neg2 = $BAMFixedData{uc($gene)}{strand}{Neg}; $neg2 = 0 if not defined $neg2;
	my $posneg = $BAMFixedData{uc($gene)}{change}{Pos}{Neg}; $posneg = 0 if not defined $posneg;
	my $negpos = $BAMFixedData{uc($gene)}{change}{Neg}{Pos}; $negpos = 0 if not defined $negpos;
	my $pospos = $BAMFixedData{uc($gene)}{change}{Pos}{Pos}; $pospos = 0 if not defined $pospos;
	my $negneg = $BAMFixedData{uc($gene)}{change}{Neg}{Neg}; $negneg = 0 if not defined $negneg;
	my $posorig = "$origDir/$gene\_Pos.orig"; ($posorig) = -e $posorig ? `wc -l $posorig` =~ /^\s*(\d+)/ : 0;
	my $negorig = "$origDir/$gene\_Neg.orig"; ($negorig) = -e $negorig ? `wc -l $negorig` =~ /^\s*(\d+)/ : 0;
	my $unkorig = "$origDir/$gene\_Unk.orig"; ($unkorig) = -e $unkorig ? `wc -l $unkorig` =~ /^\s*(\d+)/ : 0;
	$dataorig{uc($gene)}{posorig} = $posorig;
	$dataorig{uc($gene)}{negorig} = $negorig;
	$dataorig{uc($gene)}{unkorig} = $unkorig;
	my $include_cpg = $BAMData{uc($gene)}{include_cpg};
	my @print0curr = ($gene,$total,$pos,$neg,$pos2,$neg2,$pospos,$negneg,$posneg,$negpos,$posorig,$negorig,$unkorig,$include_cpg);
	for (my $i = 0; $i < @print0curr; $i++) {
		$print0curr[$i] = "NA" if not defined $print0curr[$i]; 
		$print0curr[$i] = "NA" if $print0curr[$i] =~ /^\s*$/;
	}
	$print0{gene}{uc($gene)} = join("\t", @print0curr);
	
}


sub parse_BAMFile {
	my ($BAMFile, $outLog) = @_;
	my %data;
	my $in;
	if ($BAMFile =~ /.gz$/) {
		open ($in, "zcat < $BAMFile|") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	elsif ($BAMFile =~ /.sam$/) {
		open ($in, "<", "$BAMFile") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	else {
		open ($in, "samtools view $BAMFile|") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	my $include_cpg = $BAMFile =~ /PCB\d+_bcBC\d+_plasmid/ ? "TRUE" : "FALSE";
	while (my $line = <$in>) {
		chomp($line);
		my @arr = split("\t", $line);
		next if @arr < 6;
		my ($read, $flag, $gene, $pos, $qual) = @arr;
		my ($flag2, $type);
		if ($arr[1] !~ /^(0|16)$/) {
			($read, $type, $flag, $flag2, $gene) = @arr;
		}
		my $strand = $flag eq 16 ? "Neg" : $flag eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$BAMFile$N, flag=$flag, line:\n\n$line\n\n");
		$data{uc($gene)}{$strand} ++;
		$data{uc($gene)}{total} ++;
		$data{uc($gene)}{include_cpg} = $include_cpg;
		if ($arr[1] !~ /^(0|16)$/) {
			my $strand2 = $flag2 eq 16 ? "Neg" : $flag2 eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$BAMFile$N, flag2=$flag2, line:\n\n$line\n\n");
			$data{uc($gene)}{strand}{$strand2} ++;
			$data{uc($gene)}{change}{$strand}{$strand2} ++;
		}
	}
	foreach my $gene (sort keys %data) {
		foreach my $strand2 (sort keys %{$data{uc($gene)}{strand}}) {
			LOG($outLog, date() . "${LPR}gene info$N: gene=$gene strand=$strand2 total=$data{uc($gene)}{strand}{$strand2}\n","NA");
		}
	}
	close $in;
	return \%data;
}

# get fixed sam file
my ($fixedBAMFile) = $opts->{footLoop}{BAMFixed};



my (@TXTFile) = <$footPeakFolder/99_FOOTSTATS/.PEAKSTATSTEMP/.0_*.TXT>;
if (@TXTFile == 0) {@TXTFile = <$footPeakFolder/.0_*.TXT>;}

my $data;
my $TXTFileInd = 0;
my $totalTXTFile = @TXTFile;
foreach my $TXTFile (@TXTFile[0..@TXTFile-1]) {
	my ($folder1, $fileName1) = getFilename($TXTFile, "folderfull");
	my ($footName) = $fileName1 =~ /.0_RESULTS_(.+).TXT$/;
	if ($TXTFileInd <= 20 or ($TXTFileInd > 20 and $TXTFileInd % 100 == 0)) {
		LOG($outLog, date() . "${LPR}Parsing TXTFile$N $LGN$TXTFileInd/$totalTXTFile$N $LCY$fileName1$N footName = $LCY$footName$N\n");
	}
	else {
		LOG($outLog, date() . "${LPR}Parsing TXTFile$N $LGN$TXTFileInd/$totalTXTFile$N $LCY$fileName1$N footName = $LCY$footName$N\n","NA");
	}
	my $parseName = parseName($footName . "_CH");
	my @header = qw(label gene strand window thres conv pcb bc plasmid desc);
	my ($label, $gene, $strand, $window, $thres, $conv, $pcb, $bc, $plasmid, $desc) = @{$parseName->{array}};
	$conv = "";
	$data = parseTXT($data, $TXTFile, $label, $gene, $strand, $window, $thres);
	$TXTFileInd ++;
}

my %print1;
my %print;
my $LABEL;

my $outPeakStatsFile = "$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT";
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "\n${LPR}Printing each gene into $LCY$outPeakStatsFile$N\n");
open (my $out, ">", $outPeakStatsFile) or DIELOG($outLog, date() . " Failed to write to $outPeakStatsFile: $!\n");
print $out "type\tlabel\tgene\tread_strand\twindow\tthres\tconv\tread_unique_total\tread_unique_peak_total\tread_unique_peak_fraction\tfootPeakFolder\tpeakFile\tinclude_cpg\n";
my $outPeakStatsInd = 0;
foreach my $label (sort keys %{$data}) {
	foreach my $gene (sort keys %{$data->{$label}}) {
		foreach my $read_strand (sort keys %{$data->{$label}{uc($gene)}}) {
			foreach my $window (sort keys %{$data->{$label}{uc($gene)}{$read_strand}}) {
				foreach my $thres (sort keys %{$data->{$label}{uc($gene)}{$read_strand}{$window}}) {
					foreach my $conv (sort keys %{$data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}}) {
						my $strand = $coor->{uc($gene)}{strand}; $strand = "Pos" if $strand eq "+"; $strand = "Neg" if $strand eq "-";
						my $goodconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "CG" : $read_strand eq "Neg" ? "GC" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						$goodconv = $goodconv ne "H" ? $goodconv : $read_strand eq "Pos" ? "CH" : $read_strand eq "Neg" ? "GH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						my $goodconv2 = $goodconv eq "CH" ? "CG" : $goodconv eq "GH" ? "GC" : $goodconv eq "CG" ? "CH" : "GH";
						my $revconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "GC" : $read_strand eq "Neg" ? "CG" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						$revconv = $revconv ne "H" ? $revconv : $read_strand eq "Pos" ? "GH" : $read_strand eq "Neg" ? "CH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						my $type;
						if ($outPeakStatsInd <= 100 or ($outPeakStatsInd > 100 and $outPeakStatsInd % 1000 == 0)) {
							LOG($outLog, date() . "Done $LGN$outPeakStatsInd$N: label=$label gene=$gene read_strand = $read_strand gene_strand = $strand goodconv = $goodconv\n");
						}
						else {
							LOG($outLog, date() . "Done $LGN$outPeakStatsInd$N: label=$label gene=$gene read_strand = $read_strand gene_strand = $strand goodconv = $goodconv\n","NA");
						}
						$outPeakStatsInd ++;

						my $include_cpg = $label =~ /^PCB\d+$/ ? "FALSE" : "TRUE";
						#my $cpg2 = $label =~ /^PCB\d+$/ ? "_C" : "";
						$type = 	$read_strand eq "UNK" ? "UNK" :
									($strand eq $read_strand and $conv eq $goodconv) ? "PEAK" :
									($strand ne $read_strand and $conv eq $goodconv) ? "PEAK_TEMP" : 
									($strand eq $read_strand and $conv eq $revconv) ? "PEAK_RCONV" :
									($strand ne $read_strand and $conv eq $revconv) ? "PEAK_TEMP_RCONV" :
									($strand eq $read_strand and $conv eq $goodconv2) ? "PEAK_OTHERC" :
									($strand ne $read_strand and $conv eq $goodconv2) ? "PEAK_TEMP_OTHERC" : "OTHERS";
						my ($currfootPeakFolder) = getFullpath($data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{footPeakFolder});
						my ($currpeakFile) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{peakFile};
						my ($read_unique_total) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_total};
						my ($read_unique_peak_total) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_total};
						my ($read_unique_peak_fraction) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_fraction} / 100;
						$print1{uc($gene)}{$type} += $read_unique_peak_total;
						$print1{uc($gene)}{strand} = $strand;
						$print1{uc($gene)}{OTHERS_TOTAL} ++ if $type eq "OTHERS";
						
						$print{$type} .= "$type\t$label\t$gene\t$read_strand\t$window\t$thres\t$conv\t$read_unique_total\t$read_unique_peak_total\t$read_unique_peak_fraction\t$currfootPeakFolder\t$currpeakFile\t$include_cpg\n";
						$print{ALL} = "$type\t$label\t$gene\t$read_strand\t$window\t$thres\t$conv\t0\t0\t0\t$currfootPeakFolder\t$currpeakFile\t$include_cpg\n";
						$LABEL = $label;
					}
				}
			}
		}
	}
}

open (my $out2, ">", "$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT") or die;
print $out2 "label\tstrand\t$print0{header}";
my @types = qw(PEAK PEAK_TEMP PEAK_RCONV PEAK_TEMP_RCONV OTHERS); 
for (my $i = 0; $i < @types; $i++) { 
	if (not defined $print{$types[$i]}) {
		$print{$types[$i]} = $print{ALL};
	}
	print $out $print{$types[$i]};
	print $out2 "\t$types[$i]\t$types[$i]/total";
}
print $out2 "\n";

foreach my $gene (sort keys %{$print0{gene}}) {
	my $print0 = $print0{gene}{uc($gene)};
	my $strand2 = $coor->{uc($gene)}{strand}; 
	my $strand = $strand2 eq "+" ? "Pos" : $strand2 eq "-" ? "Neg" : "UNKSTRAND";
	die "Undefined strand (strand=$strand, strand2=$strand2) at gne=$gene\n" if $strand eq "UNKSTRAND";
	my $posorig = $dataorig{uc($gene)}{posorig};
	my $negorig = $dataorig{uc($gene)}{negorig};
	my $unkorig = $dataorig{uc($gene)}{unkorig};
	print $out2 "$LABEL\t$strand\t$print0";
	for (my $i = 0; $i < @types-1; $i++) { 
		my $print1 = $print1{uc($gene)}{$types[$i]}; $print1 = 0 if not defined $print1;
		my $print1denomPos = $strand eq "Pos" ? $posorig : $strand eq "Neg" ? $negorig : $unkorig; $print1denomPos = 1 if $print1denomPos == 0;
		my $print1denomNeg = $strand eq "Pos" ? $negorig : $strand eq "Neg" ? $posorig : $unkorig; $print1denomNeg = 1 if $print1denomNeg == 0;
	#	print "type=$types[$i] strand=$strand print1=$print1 denom = $print1denomPos neg = $print1denomNeg\n";
		print $out2 "\t$print1";
		if ($i < @types - 1) {
			print $out2 "\t" . int(100 * $print1 / $print1denomPos +0.5)/100 if $i % 2 == 0;
			print $out2 "\t" . int(100 * $print1 / $print1denomNeg +0.5)/100 if $i % 2 == 1;
			LOG($outLog, "strand=$strand, $print1denomPos\n","NA") if $i % 2 == 0;
			LOG($outLog, "strand=$strand, $print1denomNeg\n","NA") if $i % 2 == 1;
		}
	}
	my $printOTHERS = $print1{uc($gene)}{OTHERS}; $printOTHERS = 0 if not defined $printOTHERS;
	my $printOTHERStot = $print1{uc($gene)}{OTHERS_TOTAL}; $printOTHERStot = 1 if not defined $printOTHERStot;
	
	my $printOTHERSperc = int($printOTHERS / $printOTHERStot * 10)/10;
	print $out2 "\t$printOTHERS\t$printOTHERSperc\n";
}
LOG($outLog, date() . "$LCY$opt_n$N Done!\n\nOutput:
$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT
$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT\n");


sub parseTXT {
	my ($data, $TXTFile, $label, $gene, $strand, $window, $thres) = @_;
	my @header = qw(footPeakFolder peakFile gene conv read_unique_total read_unique_peak_total read_unique_peak_fraction peakType);

	open (my $in1, "<", $TXTFile) or die "Cannot read from $TXTFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		if ($line =~ /^#/) {
			next;
		}
		my @arr = split("\t", $line);
		if (@header != @arr) {
			LOG($outLog, "total number of header isn't the same as total column in file $LCY$TXTFile$N!\n\n$LPR$line$N\n\n");
			my $max = @header > @arr ? @header : @arr;
			for (my $i = 0; $i < $max; $i++) {
				my $headerVal = defined $header[$i] ? $header[$i] : "UNDEF";
				my $arrVal = defined $arr[$i] ? $arr[$i] : "UNDEF";
				print "$i header=$headerVal val=$arrVal\n";
			}
			DIELOG($outLog, "\n");
		}
#		print "\n$line\n";
		my $conv = $arr[3];
		for (my $i = 0; $i < @arr; $i++) {
			next if $i == 3 or $i == 2;
			$data->{$label}{uc($gene)}{$strand}{$window}{$thres}{$conv}{$header[$i]} = $arr[$i];
		}
	}
	close $in1;
	return($data);
}

=cut
__END__
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
	  $coor{uc($gene)}{chr} = $chr;
	  $coor{uc($gene)}{beg} = $beg;
	  $coor{uc($gene)}{end} = $end;
	  $coor{uc($gene)}{strand} = $strand;
	  #print "chr=$chr:$beg-$end gene=$gene, strand=$strand\n";
   }
   close $in;
   return \%coor;
}

__END__
my %data;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

170802_pcb05_RnaseHneg_ccs_3minFP.fastq.gz.rmdup.fq.gz_PEAK_20W35T	  CALM3_Neg_20_0.35_GC.PEAK	   CALM3   GC	  160	 3	   1.9
