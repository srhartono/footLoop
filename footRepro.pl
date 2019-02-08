#!/usr/bin/perl
# 5.1

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_N $opt_1 $opt_2 $opt_m $opt_i $opt_c $opt_n $opt_G $opt_v $opt_o); #$opt_G $opt_v $opt_n $opt_r $opt_R $opt_B $opt_c); #v $opt_x $opt_R $opt_c $opt_t $opt_n);
getopts("cn:G:vo:N:1:2:i:m:");#vg:w:G:r:R:B:c");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
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

die "\nUsage: $YW$0$N -n$LCY <footPeak output directory>$N $LPR <original .fq files>$N

${LGN}Optionals$N:

-m: threshold for min read in a cluster [30]
with readID, 1, and 2 files:
-1: PCBs to compare (group 1)
-2: PCBs to compare (group 2)
-G $LGN<gene to process>$N] 
-c: This is CpG (in vitro)
-o: output

" if not defined $opt_n;
my ($min_read_in_cluster) = defined $opt_m ? $opt_m : 30;
die "\nERROR: -n footPeak dir $LCY$opt_n$N doesn't exists!\n\nUsage: $YW$0$N -n <footPeak output directory>\n\n" if not -d $opt_n;
$opt_o = $opt_n if not defined $opt_o;
#die "\nERROR: -o output dir not defined!\n" if not defined $opt_o;

my $resDir = $opt_o . "/FOOTREPRO/";
makedir($resDir);
#my $outLogFile = "$footPeakFolder/footPeak_repro_logFile.txt";
my $outLogFile = "$resDir/footPeak_repro_logFile.txt";
my ($currMainFolder) = `pwd`; chomp($currMainFolder);
my ($user) = $homedir =~ /home\/(\w+)/;
$user = "USER" if not defined $user;
my $uuid = getuuid();
my $date = date();
my $footPeakFolder = $opt_n;
my $genewant = $opt_G if defined $opt_G;

my ($OUTDIRS) = makeOutDir($resDir);
open (my $outLog, ">", $outLogFile) or die "Failed to write to $outLogFile: $!\n";
my $COOR; 

### Make sure group 1 and 2 exists ###
my %read;
my (@grp1, @grp2);
if (defined $opt_1 or defined $opt_2 or defined $opt_i) {
	if (not defined $opt_1 or not defined $opt_2) {
		my $opt1 = defined $opt_1 ? $opt_1 : "UNDEF";
		my $opt2 = defined $opt_2 ? $opt_2 : "UNDEF";
		DIELOG($outLog, date() . " Error: group 1 exists (-1 $opt1) but group 2 (-2 $opt2) doesn't!\n") if not defined $opt_2;
		DIELOG($outLog, date() . " Error: group 2 exists (-2 $opt1) but group 1 (-1 $opt1) doesn't!\n") if not defined $opt_1;
	}
	@grp1 = split(",", $opt_1);
	@grp2 = split(",", $opt_2);
}
=comment
my ($readIDFile) = "$opt_n/.readID";
$readIDFile = $opt_i if defined $opt_i;
LOG($outLog, "\n\n-----------------\n${LPR}-1. Getting group1 (and group2) from $readIDFile$N\n");
DIELOG($outLog, date() . "Cannot find readIDFile ($readIDFile)\n") if not -e $readIDFile;
open (my $inID, $readIDFile) or DIELOG($outLog, "Cannot read from $readIDFile: $!\n");
while (my $line = <$inID>) {
	chomp($line);
	$line =~ s/\t#.+//;
	my ($readID, $num, $desc) = split("\t", $line);
	$desc = "" if not defined $desc;
	$read{id}{$num} = $readID;
	$read{desc}{$num} = (not defined $desc) ? "" : $desc =~ /^#/ ? "" : $desc;
	if (not defined $opt_1) {
		push(@grp1, $num) if not grep(/^$num$/, @grp1);
	}
	LOG($outLog, "Parsed $readID\t$num\t$desc\n","NA");
}

close $inID;
if (defined $opt_1 or @grp1 != 0) {
	LOG($outLog, "Group 1:\n") if defined $opt_1;
	foreach my $num (@grp1[0..@grp1-1]) {
		my $number = $read{id}{$num};
		my $desc   = $read{desc}{$num};
		DIELOG($outLog, "\tGroup 1: Cannot find readID for number=$num\n") if not defined $number;
		$read{1}{$number}{desc} = $desc;
		$read{1}{$number}{num} = $num;
		my $descprint = $desc eq "" ? " (No Description)" : " description=$desc";
		LOG($outLog, "\t$num=$number$descprint\n");
	}
}
if (defined $opt_2) {
	LOG($outLog, "Group 2:\n");
	foreach my $num (@grp2[0..@grp2-1]) {
		my $number = $read{id}{$num};
		my $desc   = $read{desc}{$num};
		DIELOG($outLog, "\tGroup 2: Cannot find readID for number=$num\n") if not defined $number;
		$read{2}{$number}{desc} = $desc;
		$read{2}{$number}{num} = $num;
		my $descprint = $desc eq "" ? " (No Description)" : " description=$desc";
		LOG($outLog, "\t$num=$number$descprint\n");
	}
}
=cut
###0. GET READ NAMES
=comment
my %read;
if (@ARGV != 0) {
	my @fastqFolder = @ARGV;
	foreach my $fastqFolder (@ARGV) {
		my @fastq = (<$fastqFolder/*.f*q>,<$fastqFolder/*.f*q.gz>);
		LOG($outLog, date() . "Parsing fastq files from $LPR$fastqFolder$N\n");
		foreach my $fq (@fastq) {
			my ($fqName) = getFilename($fq);
			LOG($outLog, date() . "      - $LCY$fqName$N\n");
			DIELOG($outLog, date() . " Fastq $LCY$fq$N does not exist!\n") if not -e $fq;
			open (my $in, "zcat $fq|") or DIELOG($outLog, date() . "Failed to read from $fq: $!\n");
			my $linecount = 0;
			while (my $line = <$in>) {
				chomp($line);
				next if $line !~ /^\@m\d+/;
				my ($num1, $num2, $num3, $num4) = $line =~ /^.*\@m?(\d+_\d+)_(.+)?\/(\d+)\/(ccs|\d+_\d+|.*)/;
				if (not defined $num1) {
				   ($num1, $num3) = $line =~ /^.*\@m?(\d+_\d+)\/(\d+)\//;
					DIELOG($outLog, "Undefined mNNN_NNN and ZMW number at: $line\n") if not defined $num1;
					$num2 = "NA";
					$num4 = "ccs";
				}
				$num1 =~ s/_//g;
				$num4 = $num4 eq "ccs" ? 0 : $num4; $num4 =~ s/_//g;
				die "Undefined num2 at line = $line\n" if not defined $num3;
				my $number = "$num1$num3$num4";
				DIELOG($outLog, "fqfile name parse error at line=$line\n") if not defined $number;
				$read{$number}{num1} = "$num1";
				$read{$number}{num2} = "$num3";
				($fqName) = $fqName =~ /.*PCB\d+_(bcBC\d+)_/ if $fqName =~ /.*PCB\d+_bcBC\d+_/;
				($fqName) = $fqName =~ /.*(PCB\d+)/ if $fqName =~ /.*PCB\d+/ and $fqName !~ /.*PCB\d+_bcBC\d+_/;
				$read{$number}{fq} = "$fqName";
				print "\t$num1 $num3 $number\n" if $linecount < 2;
				$linecount ++;
			}
			close $in;
		}
	}
}
=cut

###0. GETTING PCB IDS FROM ONLINE DATASET
LOG($outLog, "\n\n-----------------\n${LPR}0. Getting pcb ids and readname from amazon s3 pcb_readname.tsv$N\n");
my $pcbs = myFootLib::get_pcb_readname($outLog);
DIELOG($outLog, "Failed to get pcb ids and readname\n") if not defined $pcbs;
my %pcb = %{$pcbs}; my $currcount = 0;
foreach my $number (sort keys %{$pcb{readName}}) {
	$currcount ++;
	my $PCBID = $pcb{readName}{$number};
	my ($ID) = $PCBID =~ /^PCB(\d+)$/; die if not defined $ID;
	if (not defined $opt_1) {
		LOG($outLog, "$currcount. Parsed $number: $pcb{readName}{$number}\n");
		$read{1}{$number}{pcbid} = $PCBID;
		LOG($outLog, "\t$PCBID, ID=$ID, number=$number\n");
	}
	else {
		LOG($outLog, "$currcount. Parsed $number: $pcb{readName}{$number}\n");
		if (grep(/^$ID$/, @grp1)) {
			$read{1}{$number}{pcbid} = $PCBID;
			my $grppcbids = $read{1}{$number}{grppcbids};
			my @grppcbids = defined $grppcbids ? @{$grppcbids} : ();
			push(@{$read{grppcbids1}}, $PCBID) if not grep(/^$PCBID$/, @grppcbids);
			LOG($outLog, "\t$LGN group 1$N: $PCBID, ID=$ID, number=$number\n");
		}
		if (grep(/^$ID$/, @grp2)) {
			$read{2}{$number}{pcbid} = $PCBID;
			my $grppcbids = $read{2}{$number}{grppcbids};
			my @grppcbids = defined $grppcbids ? @{$grppcbids} : ();
			push(@{$read{grppcbids2}}, $PCBID) if not grep(/^$PCBID$/, @grppcbids);
			LOG($outLog, "\t$LCY group 2$N: $PCBID, ID=$ID, number=$number\n");
		}
		if (not grep(/^$ID$/, @grp2) and not grep(/^$ID$/, @grp1)) {
			LOG($outLog, "\t$GR not in any group$N: $PCBID, ID=$ID, number=$number\n");
		}
	}
}
if (defined $opt_1) {
	my $grpnames = $read{grppcbids1}; my $grpname = defined $grpnames ? join("_", @{$grpnames}) : "NOGRPNAME1";
	$grpname =~ s/_PCB/_/g;
	print "grpname1 = $grpname = @{$grpnames}\n";
	$read{grppcbid1} = $grpname;
	$grpnames = $read{grppcbids2}; $grpname = defined $grpnames ? join("_", @{$grpnames}) : "NOGRPNAME2";
	$grpname =~ s/_PCB/_/g;
	print "grpname2 = $grpname = @{$grpnames}\n";
	$read{grppcbid2} = $grpname;
}
=comment
LOG($outLog, "\nGroup 1:\n");
#foreach my $number (sort keys %{$read{1}}) {
foreach my $number (sort @grp1) {
	my $desc = $read{1}{$number}{desc};
	my $num = $read{1}{$number}{num};
	my ($num1, $num4) = $number =~ /^(.+);(.+)$/;
	my $fqName = "FQNAME";
	$num1 =~ s/_//g;
	my $pcbid = $pcb{readName}{$num1};
	$pcbid = $number if not defined $pcbid;
	#LOG($outLog, "grp1: num1 = $num1, num4=$num4, pcb=$pcbid, desc=$desc\n");
	#my $num1 = $read{1}{$number}{num1};
	#my $num2 = $read{1}{$number}{num2};
	#my $fqName = $read{1}{$number}{fq};
	LOG($outLog, date() . "Undefined pcbid for number=$number num1=$num1 num4=$num4 desc=$desc using '$number' as pcbid!\n") if not defined $pcbid;
	my $current_readID = $pcbid;# eq $fqName ? $pcbid : "$pcbid\_$fqName";
	$read{1}{$number}{pcbid} = $pcbid;
	push(@{$read{1}{$number}{grppcbids}}, $pcbid);
	$read{1}{$num1}{pcbid} = $pcbid;
	$read{1}{$number}{num4} = $num4;
	$read{1}{$num1}{num4} = $num4;
	push(@{$read{1}{$num1}{grppcbids}}, $pcbid);
	LOG($outLog, "read1 number=$number num1=$num1 pcbid=$pcbid desc=$desc num=$num\n");
}
foreach my $number (sort keys %{$read{1}}) {
	$read{1}{$number}{grppcbid} = join("_", @{$read{1}{$number}{grppcbids}});
}

if (defined $read{2}) {
	LOG($outLog, "\nGroup 2:\n");
	foreach my $number (sort keys %{$read{2}}) {
		my $desc = $read{2}{$number}{desc};
		my $num = $read{2}{$number}{num};
		my ($num1, $num4) = $number =~ /^(.+);(.+)$/;
		my $fqName = "FQNAME";
		$num1 =~ s/_//g;
		my $pcbid = $pcb{readName}{$num1};
		$pcbid = $number if not defined $pcbid;
		LOG($outLog, date() . "Undefined pcbid for number=$number num1=$num1 num4=$num4 desc=$desc using '$number' as pcbid!\n") if not defined $pcbid;
		$pcbid = $desc if not defined $pcbid;
		$read{2}{$number}{pcbid} = $pcbid;
		push(@{$read{2}{$number}{grppcbids}}, $pcbid);
		$read{2}{$number}{num4} = $num4;
		$read{2}{$num1}{num4} = $num4;
		$read{2}{$num1}{pcbid} = $pcbid;
		push(@{$read{2}{$num1}{grppcbids}}, $pcbid);
		LOG($outLog, "read1 number=$number num1=$num1 pcbid=$pcbid desc=$desc num=$num\n");
	}
	foreach my $number (sort keys %{$read{2}}) {
		$read{2}{$number}{grppcbid} = join("_", @{$read{2}{$number}{grppcbids}});
	}
}
=cut
sub det_pcb {
	my ($readID, $readhash, $grp) = @_;
	my %read = %{$readhash};
	foreach my $number (sort keys %read) {
		if ($readID =~ /^$number/) {
			my $pcbid = $read{$number}{pcbid};
			return ($number, $pcbid, 1);
		}
	}
	my ($number) = $readID =~ /^(\d{12})/;
	   ($number) = $readID if not defined $number;
	return($number, "PCB$number", 0);
}
###1. PARSING FOOTPEAK LOGFILE TO GET GENE INFO
LOG($outLog, "\n\n-----------------\n${LPR}1. Parsing footPeak logFile $LCY$footPeakFolder/${LGN}footPeak_logFile.txt$N\n");
($COOR, $outLog) = parse_footPeak_logFile("$footPeakFolder/footPeak_logFile.txt", $footPeakFolder, $opt_G, $outLog);

###2. GET FOOTCLUST BED.INDIV.CLUST FILES
LOG($outLog, "\n\n-----------------\n${LPR}2. Getting footClust bed.indiv.clust file $LCY$footPeakFolder/${LGN}FOOTCLUST/CLUST_GENOME/*.PEAK.genome.bed.indiv.clust$N\n");
my @clustFiles = <$footPeakFolder/FOOTCLUST/CLUST_GENOME/*.PEAK.genome.bed.indiv.clust>;
my $maxClustFile = @clustFiles >= 3 ? 2 : @clustFiles == 0 ? 0 : @clustFiles-1;
LOG($outLog, date() . "$LGN" . scalar(@clustFiles) . "$N footClust files found! (see logFile.txt $LCY$outLogFile$N for examples)\n");
if ($maxClustFile > 0) {
	LOG($outLog, date() . "--> First three examples:\n" . date() . "   \t --> ","NA");
	LOG($outLog, join("\n" . date() . "   \t --> ", @clustFiles[0..$maxClustFile]) . "\n","NA");
}

###3. PARSING EACH CLUST FILES AND WRITE OUTPUT TABLE
LOG($outLog, "\n\n-----------------\n${LPR}3. Parsing each bed.indiv.clust files and writing tables$N\n");
my $fileCount = 0;
my %Rscript;
for (my $f = 0; $f < @clustFiles; $f++) {
	my $clustFile = $clustFiles[$f];
	my ($clustFilename) = getFilename($clustFile);
	my $currLog = date() . "-->$YW$f$N. Parsing $LGN$clustFilename$N\n";

	# get sample flag and info
	my $parseName  = parseName($clustFile);
	my @arr = @{$parseName->{array}};
	my ($label, $gene, $readStrand, $window, $thres, $rconvType) = @arr;
	if (defined $genewant and $genewant ne $gene) {
		LOG($outLog, $currLog . date() . "\t--> Nexted as not the same as $genewant\n","NA");
		next;
	}
	my $geneStrand = $COOR->{$gene}{STRAND};
	next if defined $opt_G and $gene ne $genewant;
	my $BEG = $COOR->{$gene}{BEG};
	my $END = $COOR->{$gene}{END};
	my $LEN = $END - $BEG;
	my $flag = $readStrand eq "Unk" ? "ALL" : "PEAK";
	my $grpnamez = "ALL";
	if (defined $opt_1) {
		$grpnamez = $read{grppcbid1} . "vs" . $read{grppcbid2};
		#print "grp = $grpnamez\n";
	}
	$flag = getFlag($clustFile, $geneStrand, $readStrand, $rconvType);
	$currLog .= date() . "\t$LCY DEBUG$N: $clustFile: flag=$flag, " . join(",", @arr) . "\n";

	#if -c then take PEAK_C only. But file doesn't contain this information. So use flag.
	if ($flag ne "PEAK" and not defined $opt_c) {
		LOG($outLog, $currLog . date() . "\tNexted as flag is $flag but not CpG (invitro -c)\n","NA");
		next;
	}
	if (defined $opt_c and $flag ne "PEAK_C") {
		LOG($outLog, $currLog . date() . "\tNexted as CpG but flag is not PEAK_C (flag=$flag)\n","NA");
		next;
	}
	$fileCount ++;

	# Parse clusterFile
	LOG($outLog, date() . "${YW}3.$fileCount$N. Parsing $LGN$clustFilename$N\n");
	open (my $in, "<", $clustFile);
	my $total = 0;
	my %data;
	my @clust;
	#1610242350461581470	72794386	72794517	1	2	1	1610242350461581470.1
	my %clusts;
	my $clustallFile = $clustFile; $clustallFile =~ s/.bed.indiv.clust/.bed.clust/;
	open (my $inAll, "<", $clustallFile) or DIELOG($outLog, "Failed to read from clustallFile $clustallFile: $!\n");
	while (my $line = <$inAll>) {
		chomp($line);
		next if $line =~ /^id\t/;
		my ($chr, $beg, $end, $geneclust, $val, $strand) = split("\t", $line);
		my ($geneid, $currclust) = $geneclust =~ /^(.+)\.(\d+)$/;
		DIELOG($outLog, "Undefined geneid or currclust at geneclust = $geneclust, line=\n\n$line\n\n") if not defined $geneid or not defined $currclust;
		$clusts{$currclust}{beg} = $beg;
		$clusts{$currclust}{end} = $end;
	}
	close $inAll;
	my %datashuf;
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /^id\t/;
#		LOG($outLog, date() . "Example line:$LCY$line$N\n") if $total == 0;
		my ($name, $x, $xmax, $y, $ymax, $clust) = split("\t", $line);
		my $current_readID = $name;
		my $len = $xmax - $x;
		my ($bestPCB) = $current_readID;
		my $isgrp = 1 if defined $opt_1;
		my ($number1, $PCB1real, $good1) = det_pcb($current_readID, \%{$read{1}}, $isgrp);
		my ($number2, $PCB2real, $good2) = det_pcb($current_readID, \%{$read{2}}, $isgrp);
		my $PCB1 =  defined $opt_1 ? $read{grppcbid1} : $PCB1real;
		$PCB1 =~ s/_PCB/_/g;
		my $PCB2 = defined $opt_2 ? $read{grppcbid2} : $PCB2real;
		$PCB2 =~ s/_PCB/_/g;
		next if $good1 eq 0 and $good2 eq 0 and defined $opt_1;
		#print "bestPCB = $bestPCB number=$number1\n";
		#my $PCB1 = $number1;#$read{1}{$bestPCB}{pcbid};
		#my $PCB2 = $number2;
#		print "HERE! name=$name\n" if $name =~ /1708140058251003640/;
#		#if (not defined $read{$bestPCB}) {
#			foreach my $pcbid (sort keys %read) {
#				if ($current_readID =~ /1708140058251003640/) {
#					print "current PCBid=$current_readID, same as pcbid=$pcbid, PCB=$read{$pcbid}{pcbid}?\n";
#				}
#				if ($current_readID =~ /$pcbid/) {
#					$bestPCB = $pcbid;
#					$PCB1 = $read{$pcbid}{pcbid};
#					print "YES\n" if ($current_readID eq "1708140058251003640");# or $current_readID eq 1708140058251003640;
#					last;
#				}
##			}
#		#}
#		#else {$PCB1 = $read{$bestPCB}{pcbid};}
		#my $idtake = 12;
		#while (not defined $PCB1) {
		#	($bestPCB) = $name =~ /^(\d{$idtake})/;
		#	DIELOG($outLog, "clustFile=$LCY$clustFile$N, name=$LCY$name$N, UNDEF bestPCB\n") if not defined $bestPCB;
		#	$PCB1 = $read{$bestPCB}{pcbid};
		#	#print "bestPCB = $bestPCB\n" if $name eq 1610242350461581470 or $name eq 161024235046898060;
		#	last if $idtake < 8;
		#	last if defined $PCB1;
		#	$idtake --;
		#}
		DIELOG($outLog, "clustFile=$LCY$clustFile$N, PCB=$LCY$current_readID$N, best=$bestPCB, name=$LCY$name$N, can't find PCB for this PCBid in database at step (0) (pcb_readname.tsv from amazon s3 ${LCY}https://s3-us-west-1.amazonaws.com/muhucsc/pcb_readname.tsv$N)\n") if (not defined $PCB1);
		LOG($outLog, "current read id=$current_readID, number1=$number1, PCB1=$PCB1real, pcbgrp1=$PCB1\n") if not defined $opt_1;
		if (defined $opt_1) {
			LOG($outLog, date() . "1: Example line:$LCY$line$N\n") if $total <= 5 and $good1 ne 0;
			#LOG($outLog, "current read id=$current_readID, number1=$number1, PCB grp1=$PCB1real, pcbgrp1=$PCB1\n") if $good1 ne 0;
			LOG($outLog, date() . "2: Example line:$LCY$line$N\n") if $total <= 5 and $good2 ne 0;
			#LOG($outLog, "current read id=$current_readID, number2=$number2, PCB grp2=$PCB2real, pcbgrp2=$PCB2\n") if $good2 ne 0;
		}
		$data{$clust}{total} ++;
		if (not defined $opt_1) {
		$data{0}{$PCB1} ++;
		$data{$clust}{PCB}{$PCB1} ++;
		$total ++;
		my %clust = ("name" => $name, "id" => $PCB1, "clust" => $clust);
		push(@clust, \%clust);
		for (my $i = 0; $i < 1000; $i++) {
			my $Xrand = int(rand($LEN - $len)) + $BEG;
			my $Xrandmax = $Xrand + $len;
			my %diff; my @currclust;
			foreach my $currclust (sort keys %clusts) {
				my $currbeg = $clusts{$currclust}{beg};
				my $currend = $clusts{$currclust}{end};
				#push(@currclust, $currclust) if ($currbeg >= $Xrand and $currbeg <= $Xrandmax) or ($Xrand >= $currbeg or $Xrand <= $currend);
				next if $Xrand < $currbeg - 100 or $Xrand > $currbeg + 100 or $Xrandmax < $currend - 100 or $Xrandmax > $currend + 100;
				$diff{$currclust} = sqrt(($currbeg - $Xrand)**2 + ($currend - $Xrandmax)**2);
				#print "$Xrand-$Xrandmax vs $currbeg-$currend: $currclust: $diff{$currclust}\n";
			}
			my $best = -1; my $bestclust = -1;
			foreach my $currclust (sort {$diff{$a} <=> $diff{$b}} keys %diff) {
				$best = $diff{$currclust}; $bestclust = $currclust; last;
			}
#			if ($bestclust eq -1) {
#				$bestclust = @currclust == 0 ? -1 : $currclust[int(rand(@currclust))];
#			}
			#$bestclust = int(rand(scalar(keys %data) - 1)) + 1 if $bestclust eq -1;
			#print "best = $best ($bestclust)\n";
			$datashuf{0}{$PCB1}[$i] ++;
			$datashuf{$bestclust}{PCB}{$PCB1}[$i] ++;
			$datashuf{$bestclust}{total}[$i] ++;
		}
		}
		if (defined $opt_1 and $good1 ne 0) {
		$data{0}{$PCB1} ++;
		$data{$clust}{PCB}{$PCB1} ++;
		$total ++;
		my %clust = ("name" => $name, "id" => $PCB1, "clust" => $clust);
		push(@clust, \%clust);
		for (my $i = 0; $i < 1000; $i++) {
			my $Xrand = int(rand($LEN - $len)) + $BEG;
			my $Xrandmax = $Xrand + $len;
			my %diff; my @currclust;
			foreach my $currclust (sort keys %clusts) {
				my $currbeg = $clusts{$currclust}{beg};
				my $currend = $clusts{$currclust}{end};
				#push(@currclust, $currclust) if ($currbeg >= $Xrand and $currbeg <= $Xrandmax) or ($Xrand >= $currbeg or $Xrand <= $currend);
				next if $Xrand < $currbeg - 100 or $Xrand > $currbeg + 100 or $Xrandmax < $currend - 100 or $Xrandmax > $currend + 100;
				$diff{$currclust} = sqrt(($currbeg - $Xrand)**2 + ($currend - $Xrandmax)**2);
				#print "$Xrand-$Xrandmax vs $currbeg-$currend: $currclust: $diff{$currclust}\n";
			}
			my $best = -1; my $bestclust = -1;
			foreach my $currclust (sort {$diff{$a} <=> $diff{$b}} keys %diff) {
				$best = $diff{$currclust}; $bestclust = $currclust; last;
			}
#			if ($bestclust eq -1) {
#				$bestclust = @currclust == 0 ? -1 : $currclust[int(rand(@currclust))];
#			}
			#$bestclust = int(rand(scalar(keys %data) - 1)) + 1 if $bestclust eq -1;
			#print "best = $best ($bestclust)\n";
			$datashuf{0}{$PCB1}[$i] ++;
			$datashuf{$bestclust}{PCB}{$PCB1}[$i] ++;
			$datashuf{$bestclust}{total}[$i] ++;
		}
		}
		if (defined $opt_1 and $good2 ne 0) {
		$data{0}{$PCB2} ++;
		$data{$clust}{PCB}{$PCB2} ++;
		$total ++;
		my %clust = ("name" => $name, "id" => $PCB2, "clust" => $clust);
		push(@clust, \%clust);
		for (my $i = 0; $i < 1000; $i++) {
			my $Xrand = int(rand($LEN - $len)) + $BEG;
			my $Xrandmax = $Xrand + $len;
			my %diff; my @currclust;
			foreach my $currclust (sort keys %clusts) {
				my $currbeg = $clusts{$currclust}{beg};
				my $currend = $clusts{$currclust}{end};
				#push(@currclust, $currclust) if ($currbeg >= $Xrand and $currbeg <= $Xrandmax) or ($Xrand >= $currbeg or $Xrand <= $currend);
				next if $Xrand < $currbeg - 100 or $Xrand > $currbeg + 100 or $Xrandmax < $currend - 100 or $Xrandmax > $currend + 100;
				$diff{$currclust} = sqrt(($currbeg - $Xrand)**2 + ($currend - $Xrandmax)**2);
				#print "$Xrand-$Xrandmax vs $currbeg-$currend: $currclust: $diff{$currclust}\n";
			}
			my $best = -1; my $bestclust = -1;
			foreach my $currclust (sort {$diff{$a} <=> $diff{$b}} keys %diff) {
				$best = $diff{$currclust}; $bestclust = $currclust; last;
			}
#			if ($bestclust eq -1) {
#				$bestclust = @currclust == 0 ? -1 : $currclust[int(rand(@currclust))];
#			}
			#$bestclust = int(rand(scalar(keys %data) - 1)) + 1 if $bestclust eq -1;
			#print "best = $best ($bestclust)\n";
			$datashuf{0}{$PCB2}[$i] ++;
			$datashuf{$bestclust}{PCB}{$PCB2}[$i] ++;
			$datashuf{$bestclust}{total}[$i] ++;
		}
		}
	}
	close $in;
	LOG($outLog, date() . "Dividing cluster shufs\n");
	foreach my $clust (sort keys %datashuf) {
		next if $clust eq 0;
		foreach my $PCB1 (sort keys %{$datashuf{$clust}{PCB}}) {
			LOG($outLog, "PCB1=$PCB1\n");
			for (my $i = 0; $i < @{$datashuf{$clust}{PCB}{$PCB1}}; $i++) {
				$datashuf{0}{$PCB1}[$i] = 0 if not defined $datashuf{0}{$PCB1}[$i];
				$datashuf{$clust}{PCB}{$PCB1}[$i] = 0 if not defined $datashuf{$clust}{PCB}{$PCB1}[$i];
				die "UNDEF clust=$clust PCBid=$PCB1 i=$i datashuftotal=$datashuf{0}{$PCB1}[$i] datashufPCBid=$datashuf{$clust}{PCB}{$PCB1}[$i]\n" if not defined $datashuf{$clust}{PCB}{$PCB1}[$i] or not defined $datashuf{0}{$PCB1}[$i];
				$datashuf{$clust}{PCB}{$PCB1}[$i] = $datashuf{0}{$PCB1}[$i] == 0 ? 0 : 
int($datashuf{$clust}{PCB}{$PCB1}[$i] / $datashuf{0}{$PCB1}[$i] * 10000 + 0.5)/100;
			}
		}
	}

	# get random reads for shuffle
#	my $shufhash = get_random_reads(\%data, \@clust, $outLog);
#	DIELOG($outLog, "Failed to get shufhash from get_random_reads!\n") if not defined $shufhash;
#	my %shuf = %{$shufhash};

	#open and print to repro table, which makes corr heatmap
	LOG($outLog, date() . "PRintint output\n");
	my $forbar;
	LOG($outLog, "\n\n#clust\tcombine");
	#my $total = 0;
	my $used = 0;
	foreach my $PCB (sort keys %{$data{0}}) {
		my $total_PCBid_clust = $data{0}{$PCB};
		LOG($outLog, "PCB=$PCB total = $total_PCBid_clust\n");
		next if $total_PCBid_clust < $min_read_in_cluster;
		$used ++;
#		LOG($outLog, "\t$PCB");
	}
	LOG($outLog, "\n");
	open (my $outTableForBar, ">", "$resDir/$clustFilename.$flag.repro.$grpnamez.tsv") or DIELOG($outLog, date() . "Failed to write to $resDir/$clustFilename.$flag.repro.$grpnamez.tsv: $!\n");
	print $outTableForBar "cluster\tsample\ttotal\ttype\tmean\tse\tpval\n";
	$forbar .= "#cluster\tsample\ttype\tmean\tse\tpval\n";
	foreach my $clust (sort {$a <=> $b} keys %data) {
		next if $clust eq 0;
		my $total_clust = $data{$clust}{total};
		print $outTableForBar "$clust\tPCB0\t$total\torig\t" . int($total_clust / $total * 10000 + 0.5)/100 . "\t0\t\"\"\n";
		$forbar .= "$clust\tPCB0\t$total\torig\t" . int($total_clust / $total * 10000 + 0.5)/100 . "\t0\t\"\"\n";
	}
	foreach my $clust (sort {$a <=> $b} keys %data) {
		next if $clust eq 0;
		my $total_clust = $data{$clust}{total};
		my $pcbcount = 0;
		my $clustperc = int($total_clust / $total * 10000 + 0.5)/100;
		LOG($outLog, "$clust\t$clustperc");
		my $print0 = 0;
		foreach my $current_readID (sort keys %{$data{0}}) {
#			my $PCB = $current_readID;
			my $PCB = $current_readID;#$read{$current_readID}{pcbid};
			DIELOG($outLog, date() . "Can't find PCB of PCBid=$LGN$current_readID$N, clustFile=$LPR$clustFile$N\n\n") if not defined $PCB;
			#my $PCB = $pcb{readName}{$current_readID};
			$pcbcount ++;
			my $total_PCBid_clust = $data{0}{$current_readID};
			my $total_clust_perc = defined $data{$clust}{PCB}{$current_readID} ? int($data{$clust}{PCB}{$current_readID} / $total_PCBid_clust * 10000 + 0.5)/100 : 0;
#			my $shufmean = int(100*$shuf{$clust}{$current_readID}{mean}+0.5)/100;
#			my $shufse = int(100*$shuf{$clust}{$current_readID}{sd}+0.5)/100;
#			my $shufpval = $shuf{$clust}{$current_readID}{pval};
#			my $pvalstar = $shufpval <= 0.05 ? "*" : "\"\"";
			my $shufrandtotal = defined $datashuf{0}{$current_readID} ? int(mean(@{$datashuf{0}{$current_readID}})) : 0;
			my $shufrandmean = defined $datashuf{$clust}{PCB}{$current_readID} ? int(mean(@{$datashuf{$clust}{PCB}{$current_readID}}) * 100)/100 : 0;
			my $shufrandse = defined $datashuf{$clust}{PCB}{$current_readID} ? int(se(@{$datashuf{$clust}{PCB}{$current_readID}}) * 100)/100 : 0;
			next if $total_PCBid_clust < 30;
			LOG($outLog, "\t$total_clust_perc");
			printf $outTableForBar "$clust\t$PCB\t$total_PCBid_clust\torig\t$total_clust_perc\t0\t\"\"\n";
			printf $outTableForBar "$clust\t$PCB\t$shufrandtotal\tshufrand\t$shufrandmean\t$shufrandse\t\"\"\n";
			$forbar .= "$clust\t$PCB\t$total_PCBid_clust\torig\t$total_clust_perc\t0\t\"\"\n";
			$forbar .= "$clust\t$PCB\t$shufrandtotal\tshufrand\t$shufrandmean\t$shufrandse\t\"\"\n";
			if ($clust eq 1) {
				my $maxClust = scalar(keys %data);
#				$shufrandmean = defined $datashuf{'-1'}{PCB}{$current_readID} ? int(1/($maxClust-1) * mean(@{$datashuf{'-1'}{PCB}{$current_readID}}) * 100)/100 : 0;
#				$shufrandse = defined $datashuf{'-1'}{PCB}{$current_readID} ? int(1/($maxClust-1) * se(@{$datashuf{'-1'}{PCB}{$current_readID}}) * 100)/100 : 0;
				$shufrandmean = defined $datashuf{'-1'}{PCB}{$current_readID} ? int(mean(@{$datashuf{'-1'}{PCB}{$current_readID}}) * 100)/100 : 0;
				$shufrandse = defined $datashuf{'-1'}{PCB}{$current_readID} ? int(se(@{$datashuf{'-1'}{PCB}{$current_readID}}) * 100)/100 : 0;
				for (my $k = -1; $k > -2; $k--) {#-1 * $maxClust; $k--) {
					printf $outTableForBar "$k\tPCB0\t$total\torig\t0\t0\t\"\"\n" if $print0 == 0; $print0 = 1;
					printf $outTableForBar "$k\t$PCB\t$total_PCBid_clust\torig\t0\t0\t\"\"\n";
					printf $outTableForBar "$k\t$PCB\t$shufrandtotal\tshufrand\t$shufrandmean\t$shufrandse\t\"\"\n";
					$forbar .= "$k\tPCB0\t$total\torig\t0\t0\t\"\"\n";
					$forbar .= "$k\t$PCB\t$total_PCBid_clust\torig\t0\t0\t\"\"\n";
					$forbar .= "$k\t$PCB\t$shufrandtotal\tshufrand\t$shufrandmean\t$shufrandse\t\"\"\n";
				}
			}

#			printf $outTableForBar "$clust\t$PCB\t$total_PCBid_clust\tshuf\t$shufmean\t$shufse\t$pvalstar\n";
#			$forbar .= "$clust\t$PCB\t$total_PCBid_clust\tshuf\t$shufmean\t$shufse\t$pvalstar\n";
#			my @shufarr = @{$shuf{$clust}{$current_readID}{array}};
#			for (my $i = 0; $i < @shufarr; $i++) {
#				my $shufperc = int($shufarr[$i] * 10000 + 0.5)/100;
#				print $outTableForBar "$clust\t$PCB\t$total_PCBid_clust\tshufindiv\t$shufperc\t0\t\"\"\n";
#			}
		}
		LOG($outLog, "\n");
	}
	LOG($outLog, "\n");
	LOG($outLog, "\n$forbar\n","NA");
	close $outTableForBar;


	#####WRITE TO R SCRIPT
	my $Rfilename = "$resDir/$clustFilename.$flag.repro.cor.$grpnamez.R";
	my $Rpdfname  = "$resDir/$clustFilename.$flag.repro.cor.$grpnamez.pdf";
	my $totalclust = (keys %data) - 1;
	my $totalsample = (keys %{$data{0}});
	LOG($outLog, "total clust = $totalclust, total sample=$totalsample, used = $used\n");
	my $Rscript = ($used > 1 and $totalclust > 1 and $totalsample > 1) ? get_Rscript($resDir, $clustFilename, $totalclust, $totalsample, $flag, $grpnamez) : 
"pdf(\"$Rpdfname\");plot(NA,xlim=c(0,1),ylim=c(0,1),bty=\"n\",xlab=NA,ylab=NA,axes=F,main=\"$clustFilename\n$flag\n$totalclust clusters\n$totalsample samples\n$used uesd\");dev.off()\n";
	open (my $RscriptOut, ">", $Rfilename) or DIELOG($outLog, date() . "Failed to write to $Rfilename: $!\n");
	print $RscriptOut $Rscript;
	close $RscriptOut;

	$Rscript{$Rfilename} = $Rpdfname;
}

###4. RUN R SCRIPTS
LOG($outLog, "\n\n-----------------\n${LPR}4. Running $LGN" . scalar(keys %Rscript) . "$LPR R scripts$N\n");

foreach my $Rscript (sort keys %Rscript) {
	system("run_Rscript.pl $Rscript > $Rscript.LOG 2>&1") == 0 or DIELOG($outLog, date() . "\n\nFailed to run R script: $!\n\nrun_Rscript.pl $Rscript\n\n");
	my $RLog = `tail -n 20 $Rscript.LOG`;
	LOG($outLog, "\n\n" . date() . "last 20 lines of R Log:\n\n$RLog\n\n");
	my $pdf = getFullpath($Rscript{$Rscript});
	LOG($outLog, "${YW}scp mitochi\@crick.cse.ucdavis.edu:$pdf ./$N\n");
}

###############
# SUBROUTINES #
###############

sub get_random_reads {
	my ($data, $clusts, $outLog) = @_;
	my %data = %{$data};
	my @clust = @{$clusts};
	my %shuf;

	foreach my $current_readID (sort keys %{$data{0}}) {
		my $total_PCBid = $data{0}{$current_readID};
		for (my $shufnum = 0; $shufnum < 500; $shufnum ++) {
			@clust = shuffle(\@clust);
			my %shufcurr;
			for (my $i = 0; $i < $total_PCBid; $i++) {
				my $clust = $clust[$i]->{clust};
				$shufcurr{$clust} ++;
			}
			foreach my $clust (sort keys %data) {
				next if $clust eq 0;
				my $shufcurrtot = defined $shufcurr{$clust} ? $shufcurr{$clust} : 0;
				push(@{$shuf{$clust}{$current_readID}{array}}, $shufcurrtot / $total_PCBid);# * 100 + 0.5)/100);
				LOG($outLog, date() . "example: clust=$clust, shufcurrtotal=$shufcurrtot, total_in_clust=$total_PCBid\n", "NA") if $clust eq 1 and $shufnum == 0;
			}
		}
		foreach my $clust (sort keys %data) {
			next if $clust eq 0;
			$shuf{$clust}{$current_readID}{mean} = int(mean(@{$shuf{$clust}{$current_readID}{array}}) * 10000 + 0.5)/100;
			$shuf{$clust}{$current_readID}{sd} = int(sd(@{$shuf{$clust}{$current_readID}{array}}) * 10000 + 0.5)/100;
			my $total_PCBid_clust = $data{0}{$current_readID};
			my $total_clust_perc = defined $data{$clust}{PCB}{$current_readID} ? int($data{$clust}{PCB}{$current_readID} / $total_PCBid_clust * 10000+0.5)/100 : 0;
			my ($nge, $nle, $ntot) = (0,0,0);
			for (my $i = 0; $i < @{$shuf{$clust}{$current_readID}{array}}; $i++) {
				my $shufvalue = $shuf{$clust}{$current_readID}{array}[$i];
				$nge ++ if $shufvalue >= $total_clust_perc;
				$nle ++ if $shufvalue <= $total_clust_perc;
				$ntot ++;
			}
			my $pval = $nge <= $nle ? ($nge + 1) / ($ntot + 1) : ($nle + 1) / ($ntot + 1);
			$pval = int(10000*$pval+0.5)/10000;
			$shuf{$clust}{$current_readID}{pval} = $pval;
			LOG($outLog, date() . "DEBUG: clust=$clust, pval=$pval, mean = $shuf{$clust}{$current_readID}{mean}, totalperc = $total_clust_perc\n","NA") if $clust == 1;
		}	
	}
	return(\%shuf);
}


sub get_Rscript {
	my ($resDir, $clustFilename, $totalclust, $totalsample, $flag, $grpnamez) = @_;
	my $Rscript = "
.libPaths( c(\"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/\", \"/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/\", .libPaths()) )
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(GMD)

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

my2dec = function(x) {
   for (i in 1:dim(x)[1]) {
      for (j in 1:dim(x)[2]) {
         x[i,j] = as.integer(x[i,j] * 100 + 0.5) / 100
      }
   }
   x
}


df = read.table(\"$resDir/$clustFilename.$flag.repro.$grpnamez.tsv\",header=T,sep=\"\\t\")

df.orig = df[df\$type == \"orig\" | df\$type == \"shufrand\",]
df.orig\$sample = as.character(df.orig\$sample)
df.orig[df.orig\$type == \"shufrand\",]\$sample = paste(df.orig[df.orig\$type == \"shufrand\",]\$sample,\"_shuf\",sep=\"\")
df.orig = df.orig[order(df.orig\$cluster, df.orig\$sample),]
if(length(df.orig[df.orig\$mean == 0,]\$mean) > 0) {
	df.orig[df.orig\$mean == 0,]\$mean = 1
}
for (i in 1:length(df.orig\$mean)) {
	df.orig\$mean[i] = df.orig\$mean[i] + i / 10000000 # so stdev will never be zero..
}
df2.orig = as.data.frame(matrix(ncol=length(unique(df.orig\$sample)),nrow=length(unique(df.orig\$cluster)),df.orig\$mean,byrow=T))
colnames(df2.orig) = unique(df.orig\$sample)
rownames(df2.orig) = unique(df.orig\$cluster)
print(\"$resDir/$clustFilename.$flag.repro.cor.$grpnamez.Rdata\")
print(\"Table:\")
print(df2.orig)
save(df2.orig,file=\"$resDir/$clustFilename.$flag.repro.cor.$grpnamez.Rdata\")
print(cor(log(df2.orig)))
print(cor((df2.orig[rownames(df2.orig != -1),])))
print(cor(log(df2.orig[rownames(df2.orig != -1),])))
df2.orig = df2.orig[rownames(df2.orig != -1),]
origcor = my2dec(cor(df2.orig))
rownames(origcor) = unique(paste(df.orig\$sample,\"\\n(\",df.orig\$total,\")\",sep=\"\"))
print(\"Correlation:\")
print(origcor)
pdf(\"$resDir/$clustFilename.$flag.repro.cor.$grpnamez.pdf\",height=9,width=9);
heatmap.3(origcor,
Rowv=TRUE,
Colv=TRUE,
cluster.by.row=TRUE,
cluster.by.col=TRUE,
dendrogram=\"both\",
color.FUN=function(x)rev(brewer.pal(11,\"RdBu\")),
breaks=c(seq(-1,-0.2,by=0.2),-0.05,0.05,seq(0.2,1,by=0.2)),
main=\"$clustFilename $flag\\nCorrelation Plot (total samples=$totalsample, total clusters=$totalclust)\",
cex.main=0.7,
cexRow=0.5,
cexCol=0.5,
cellnote=origcor,
notecol=\"white\")

ggplot(df.orig,aes(as.factor(cluster),mean)) + geom_bar(aes(fill=sample),stat=\"identity\",position=\"dodge\") + theme_bw() + theme(panel.grid=element_blank()) +
	ggtitle(\"$clustFilename $flag\\nPercent Reads in each Cluster\\n(total samples=$totalsample, total clusters=$totalclust)\") +
	coord_cartesian(ylim=c(0,100)) + xlab(\"Cluster\") + ylab(\"Percent\")

write.table(origcor,file=\"$resDir/$clustFilename.$flag.repro.COR\",sep=\"\\t\",col.names=T,row.names=T,quote=F)
options(width=200)
print(origcor)
";
	return $Rscript;
}





__END__
pdf(\"Clust.pdf\",height=13*3.5,width=7)
ggplot(dm3,aes(as.factor(variable),value)) + geom_bar(aes(fill=as.factor(variable)),stat = \"identity\",position=position_dodge(width=0.75)) + theme_bw() +ylab(\"% Read\") + xlab(\"Cluster Number\") +
	geom_text(aes(group=as.factor(variable),label=star),stat=\"identity\",position=position_dodge(width=0.75)) + theme(panel.grid = element_blank()) + facet_grid(cluster~.)
dev.off()


#corrplot(cor(df2),col=rev(brewer.pal(11,\"RdBu\")),low=0,upp=1,order=\"hclust\",method=\"color\",addgrid=1,number);
#corrplot(cor(df2),col=\"black\",method=\"number\",add=T);

#pdf(\"$resDir/$clustFilename.$flag.repro.pdf\",width=7,height=15)
#temp = df[df\$type == \"orig\" | df\$type == \"shuf\",]
#p = ggplot(temp,aes(as.factor(sample),mean)) + geom_bar(aes(fill=as.factor(type)),stat = \"identity\",position=position_dodge(width=0.75)) + theme_bw() +ylab(\"% Read\") + xlab(\"Cluster Number\") +
#	geom_errorbar(aes(group=as.factor(type),ymin=mean-se,ymax=mean+se),stat=\"identity\",position=position_dodge(width=0.9),width=0.5) +
#	geom_text(aes(group=as.factor(type),label=pval),stat=\"identity\") + theme(panel.grid = element_blank()) + facet_grid(cluster~.)
#print(p)
#dev.off()

