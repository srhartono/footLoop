#!/usr/bin/perl
	
use strict; use warnings; use Getopt::Std; use File::Basename qw(dirname); use Cwd qw(abs_path);
use vars qw($opt_v $opt_b $opt_i $opt_l $opt_f);
getopts("vb:i:l:f");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite; use footPeakAddon;
use feature 'say';

my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";
my @version = `cd $footLoopDir && git log | head `;
my $version = "UNKNOWN";
foreach my $line (@version[0..@version-1]) {
   if ($line =~ /^\s+V\d+\.?\d*\w*\s*/) {
      ($version) = $line =~ /^\s+(V\d+\.?\d*\w*)\s*/;
   }
}
if (not defined $version or (defined $version and $version eq "UNKNOWN")) {
   ($version) = `cd $footLoopDir && git log | head -n 1`;
}

my ($input1) = ($opt_i);
my $barcodePCB11 = "/home/mitochi/pacbio_indexes/barcode_PCB11.tsv";
my $barcode170726 = "/home/mitochi/pacbio_indexes/barcode_170726.tsv";
my $barcode171005 = "/home/mitochi/pacbio_indexes/barcode_171005.tsv";

die "
Usage: $YW$0$N -b $LPR<barcode>$N -i $CY<folder/lbc>$N [optional: -l $LGN<label>$N]

Input:
-i not defined: process all lbc.fastq files in current folder
-i <folder>: process all lbc.fastq files in the given folder

Barcode:
-b 1: barcode_171005.tsv
-b 2: barcode_170726.tsv
-b 3: barcode_PCB11.tsv
	
" unless defined $opt_b;# $opt_i and -e $opt_i and defined $opt_b;
my ($date) = getDate();
my ($uuid) = getuuid();
my $verbose = defined $opt_v ? "" : "NA";
my @ignore = `cat IGNORE` if -e ".IGNORE";

my ($outLog, $run_script, @files) = get_input($input1);
my $data = process_barcode($files[0], $opt_b, $outLog);
my $fileCount = 0;
my %file;

foreach my $file (sort @files) {
	next if $file =~ /.rmdup.fq/;
	my ($folder1, $fileName1) = getFilename($file, "folderfull");
	if ($fileCount == 0) {
		makedir("$folder1/.originallbc/") if not -d "$folder1/.originallbc/";
		my $input = (not defined $opt_i) ? " -i " . `pwd` : " -i " . $opt_i; chomp($input);
		my $verboses = defined $opt_v ? " -v" : "";
		my $labelz = defined $opt_l ? $opt_l : "";
	}
	$fileName1 =~ s/bc[ \t,\-\_\+\=\!\@\#\$\%\^\&\*\(\)\\\/\[\]\{\}\>\<\'\"\;\:\~\`]*(\d+)/bc$1/ig;
	my ($bcNum) = $fileName1 =~ /bc(\d+)/i;
	DIELOG($outLog, "\n\nERROR: multiple BCNUM $LGN$bcNum$N exists ($file and $file{$bcNum})\n\n") if defined $file{$bcNum};
	$file{$bcNum} = $file;
	$fileCount ++;
}
LOG($outLog, "\n\n$YW ----------- Processing$LGN $fileCount$YW Fastq Files in $LCY$input1$YW ----------- $N\n\n");
$fileCount = 0;
my %sbatch;
foreach my $bcNum (sort {$a <=> $b} keys %file) {
	$fileCount ++;
	#LOG($outLog, "\n$LGN$fileCount/" . scalar(keys %file) . "$N. bcNum=$YW$bcNum$N: Processing $LCY$file{$bcNum}$N\n") if defined $opt_v;
	my ($res, $resFile) = main($file{$bcNum}, $data, $fileCount, $run_script, $outLog);
	LOG($outLog, "\n$LGN$fileCount/" . scalar(keys %file) . "$N. bcNum=$YW$bcNum$N: Processed $LCY$file{$bcNum}$N into $LPR$resFile$N\n");# if not defined $opt_v;
	LOG($outLog, "$res\n", $verbose);
	#LOG($outLog, "\n$YW$bcNum$N: Processed $LCY$file{$bcNum}$N into $LPR$res$N\n");
	#print "$res\n";
	$sbatch{$bcNum} = $resFile;
}
LOG($outLog, "\n\n$YW ----------- Running $LGN$fileCount$YW Sbatches ----------- $N\n\n");
foreach my $bcNum (sort {$a <=> $b} keys %file) {
	my $slurmid = -1;
	if (-d "$sbatch{$bcNum}.rmdup.fq.gz_MAP_Z" and -d "$sbatch{$bcNum}.rmdup.fq.gz_PEAK_t35" and not defined $opt_f) {
		LOG($outLog, "$sbatch{$bcNum}\_dedupe.sbatch #${LGN}BC$N BC$bcNum$N,${LCY}slurmID$N $slurmid ${LGN}ALREADY EXISTS$N!\n");
		next;
	}
	($slurmid) = `sbatch $sbatch{$bcNum}\_dedupe.sbatch`  =~ /job (\d+)$/; 
	if (not defined $slurmid) {
		LOG($outLog, "$sbatch{$bcNum}\_dedupe.sbatch #${LGN}BC$N BC$bcNum$N,${LCY}slurmID$N $slurmid FAILED TO RUN!\n");
		next;
	}
	chomp($slurmid); 
	LOG($outLog, "$sbatch{$bcNum}\_dedupe.sbatch #${LGN}BC$N BC$bcNum$N,${LCY}slurmID$N $slurmid\n");
}
LOG($outLog, "\nUUID $uuid Done\n");



###############
# SUBROUTINES #
###############

sub get_input {
	my ($input1) = @_;
	$input1 = "./" if (not defined $input1);
	print "\n\nERROR! file $LCY$input1$N does not exists!\n\n" and die if not -e $input1;
	my $inputFullpath = getFullpath($input1);
	my ($folder1, $fileName1) = getFilename($input1, "folderfull");
	my ($folder2, $fileName2) = getFilename($inputFullpath, "folderfull");
	$folder2 = $inputFullpath if -d $inputFullpath;
	$folder1 = "./" if $folder1 eq "";
	makedir("$folder1/originallbc/") if not -d "$folder1/.originallbc/";
	my $input = (not defined $opt_i) ? " -i " . `pwd` : " -i " . $opt_i; chomp($input);
	my $verboses = defined $opt_v ? " -v" : "";
	my $labelz = defined $opt_l ? " -l $opt_l" : "";

	my $run_script = "$0 -b $opt_b$verboses$labelz -i $inputFullpath";
	open ($outLog, ">", "$folder1/lbc2barcode_logFile.txt") or DIELOG($outLog, "\n\nFailed to write to $folder1/lbc2barcode_logFile.txt: $!\n\n");
	LOG($outLog, ">$0 Version $version\n>UUID: $uuid\n>Date: $date\n>Run script: $0 -b $opt_b$verboses$labelz -i $inputFullpath\n>Folder: $folder2\n");
	LOG($outLog, "\n\n$YW ----------- Getting input files from $LCY$input1$YW ----------- $N\n\n");
	my (@files) = -d $input1 ? (<$input1/*lbc*.f*q>, <$input1/*lbc*.f*q.gz>, <$input1/*BC*.f*q>, <$input1/*BC*.f*q.gz>) : ($input1);
	@files = sort @files;
	DIELOG($outLog, "\n\nError! There is zero lbc*.f*q files in folder $LCY$input1$N\n\n") if @files == 0;
	LOG($outLog, date() . "Input Files:\n");
	for (my $i = 0; $i < @files; $i++) {
		LOG($outLog, date() . "$LGN$i$N. $LCY$files[$i]$N\n");
	}
	return ($outLog, $run_script, @files);
}

sub process_barcode {
	my ($input1, $optb, $outLog) = @_;
	my $PCB = $opt_l;
	if (not defined $opt_l or $opt_l eq "") {
		($PCB) = getFullpath($input1);
	}
	my $PCB0 = $PCB;
	$PCB =~ s/^.*PCB[\t \-\,\.\_\=\+]*(\d+).*$/PCB$1/i;
	DIELOG($outLog, "1. Cannot parse PCB from $PCB0 (make sure to use -l PCB followed by digits, e.g. PCB15)\n") if not defined $PCB or (defined $PCB and $PCB !~ /^PCB\d+$/);
	$PCB = uc($PCB);
	my $barcodeFile = $optb eq 1 ? $barcode171005 : $optb eq 2 ? $barcode170726 : $optb eq 3 ? $barcodePCB11 : $optb;
	print "\n\nERROR! barcodeFile $LCY$barcodeFile$N does not exist!\n\n" and die if not -e $barcodeFile;
	my %data;
	open (my $in1, "<", $barcodeFile) or die "Cannot read from $barcodeFile: $!\n";
	LOG($outLog, "\n\n$YW ----------- Parsing Barcode $LCY$barcodeFile$YW ----------- $N\n\n");
	LOG($outLog, date() . "Barcodes:\n");
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		next if $line =~ /^No/ or $line =~ /Barcode/i;
		my ($no, $junk, $len, $conc, $vol, $label2, $barcode, $desc, $plasmid, $primer1, $primer2) = split("\t", $line);
		die "Failed to parse desc from line\n\n$line\n\n" if not defined $primer2 or not defined $desc;
		$desc = uc($desc);
		$desc =~ s/[ \t,\-\_\+\=\!\@\#\$\%\^\&\*\(\)\\\/\[\]\{\}\>\<\'\"\;\:\~\`]+//g;
		die "Failed to parse plasmid from line\n\n$line\n\n" if not defined $primer2 or not defined $plasmid;
		$plasmid = uc($plasmid);
		$plasmid =~ s/[ \t,\-\_\+\=\!\@\#\$\%\^\&\*\(\)\\\/\[\]\{\}\>\<\'\"\;\:\~\`]+//g;
		my ($bcNum) = $label2 =~ /^BC(\d+)$/;
		die "Failed to parse bcNum number from label2=$LGN$label2$N\n\n" if not defined $bcNum;
		my $lbc = "lbc$bcNum--lbc$bcNum";
		my $labelz = "$PCB\_$label2\_$plasmid\_$desc";
		$data{$lbc}{label} = "$labelz";
		$data{$lbc}{barcode} = $barcode;
		$data{$lbc}{plasmid} = $plasmid =~ /^PFC(8|19)$/i ? "pFC8_fixed.bed" : $plasmid =~ /^PFC53$/i ? "pFC53_fixed.bed" : $plasmid =~ /^PFC66$/i ? "pFC66_fixed.bed" : $plasmid =~ /^PBSPHIX$/i ? "PBSPHIX.bed" : "PLASMIDUNK";
		$data{$lbc}{dedupe} = "dedupe.sh ow=t k=31 nam=4 mo=200 e=30 mid=98 in=\"DEDUPEINPUT\" out=\"DEDUPEOUTPUT\" && gzip -f DEDUPEOUTPUT && md5sum DEDUPEOUTPUT.gz > DEDUPEOUTMD5\n";
		$data{$lbc}{map} = "footLoop.pl -r DEDUPEOUTPUT.gz -Z -n FOOTLOOPOUTPUT_MAP_Z -i /home/mitochi/pacbio_indexes/$data{$lbc}{plasmid} -g /home/mitochi/pacbio_indexes/plasmid_fixed.fa -l $labelz\n";
		$data{$lbc}{peak} = "footPeak.pl -n FOOTLOOPOUTPUT_MAP_Z -o FOOTPEAKOUTPUT_PEAK_t35 -t 35 -l $labelz\n";
		$data{$lbc}{peak} .= "footPeak.pl -n FOOTLOOPOUTPUT_MAP_Z -o FOOTPEAKOUTPUT_PEAK_t55 -t 55 -l $labelz\n";
		$data{$labelz}{label} = $labelz;
		$data{$labelz}{barcode} = $data{$lbc}{barcode};
		$data{$labelz}{plasmid} = $data{$lbc}{plasmid};
		$data{$labelz}{dedupe} = $data{$lbc}{dedupe};
		$data{$labelz}{map} = $data{$lbc}{map};
		$data{$labelz}{peak} = $data{$lbc}{peak};
		LOG($outLog, date() . "$LGN$no$N. bcNum=$LGN$bcNum$N, lbc=$LPR$lbc$N, label=$LGN$labelz$N\n");
	}
	close $in1;
	return \%data;
}	

sub main {
	my ($input, $data, $num, $run_script, $outLog) = @_;


	my $verbose = "NA" if not defined $opt_v;
	my ($input1) = getFullpath($input);
	my ($label) = $input1 =~ /(PCB[\,\.\-\_]*\d+)/i;
	$label = $opt_l if defined $opt_l;
	DIELOG($outLog, "Cannot parse label from $input1!\n") if not defined $label;
	my ($folder1, $fileName1) = getFilename($input, "folderfull");
	$fileName1 =~ s/BC[ \t,\-\_\+\=\!\@\#\$\%\^\&\*\(\)\\\/\[\]\{\}\>\<\'\"\;\:\~\`]*(\d+)/BC$1/g;
	my ($bcNum) = $fileName1 =~ /bc(\d+)/i;
	my ($lbc, $ext);
	($lbc, $ext) = $fileName1 =~ /^(lbc\d+\-\-lbc\d+)\.(f.*q.gz|f.*q|f.*q.gz.*)$/i;
	($lbc, $ext) = $fileName1 =~ /^($label\_BC\w+)\.(f.*q.gz|f.*q|f.*q.gz.*)$/i if not defined $lbc or not defined $ext;
	DIELOG($outLog, "\n\nERROR: input1 $LGN$input1$N with fileName1 $LCY$fileName1$N does not contain lbc or ext in its name!\n") and die if not defined $ext or not defined $lbc or not defined $bcNum;
#	LOG($outLog, "\n\nERROR: input1 $LGN$input1$N with fileName1 $LCY$fileName1$N does not contain lbc in its name!\n") and die if not defined $bcNum;
	$bcNum = uc($bcNum);
	$ext =~ s/fastq/fq/;
	

#	foreach my $lbc (sort keys %{$data}) {
#		next if $lbc =~ /^BC/;
	my $res = "";
	#$ext = "" if not defined $ext;
	#$fileName1 =~ s/lbc[ \t,\-\_\+\=\!\@\#\$\%\^\&\*\(\)\\\/\[\]\{\}\>\<\'\"\;\:\~\`]*(\d+)/lbc$1/g;
	#my $lbc = $fileName1;#"lbc$bcNum--lbc$bcNum";
	LOG($outLog, date() . "$LGN$num$N. bcNum=$LGN$bcNum$N, filename=$LCY$fileName1$N, lbc=$lbc, ext=$ext\n", $verbose);
	if (not defined $data->{$lbc}) {
		$res .= "\n\nERROR: SKIPPED file $LRD$input1$N cannot find its barcode num $LGN$bcNum$N ($LGN$lbc$N) in barcode file -b $YW$opt_b$N\n\n";
		return($res, "(FILE SKIPPED)");
	}
	my $labelz = $data->{$lbc}{label};
	my $barcode = $data->{$lbc}{barcode};
	my $resFile = "$folder1/$labelz.$ext";
	my $resrmdupMD5  = "$folder1/.$labelz.$ext.rmdup.fq.gz.md5";
	$data->{$lbc}{lbcpreproc} = "";
	my ($folder0, $fileName0) = getFilename($input1, "folderfull");
	if ($lbc eq "lbc$bcNum--lbc$bcNum") {
		$res .= date() . "\tbcNum=$LGN$bcNum$N, lbc=$LPR$lbc$N, label=$LGN$labelz$N\n";
		$data->{$lbc}{lbcpreproc} = "/bin/cp $input1 $folder1/.originallbc/\n";
		$res .= date() . "\t/bin/cp $LCY$input1$N $LGN$folder1/.originallbc/$N\n";#$verbose);
		$data->{$lbc}{lbcpreproc} .= "/bin/mv $input1 $folder1/$labelz.$ext\n";
		$res .= date() . "\t/bin/mv $LCY$input1$N $LGN$folder1/$labelz.$ext$N\n";
		if ($ext !~ /gz$/) {
			my $origFile = "$folder1/.originallbc/$fileName0";
			my $origFileGz = "$folder1/.originallbc/$fileName0.gz";
			my $origFileGzMd5 = "$folder1/.originallbc/.$fileName0.gz.md5";
			#gzip
			$res .= date() . "\tgzip -f $LCY$origFile$N\n";
			#md5
			$res .= date() . "\tmd5sum $LCY$origFileGz$N > $YW$origFileGzMd5$N\n";
			$data->{$lbc}{lbcpreproc} .= "gzip -f $origFile\n";
			$data->{$lbc}{lbcpreproc} .= "md5sum $origFileGz > $origFileGzMd5\n";
		}
	}
	if ($ext !~ /gz$/) {
		$resFile .= ".gz";
		$resrmdupMD5  = "$folder1/.$labelz.$ext.gz.rmdup.fq.gz.md5";
		$res .= date() . "\tgzip -f $LCY$folder1/$labelz.$ext$N\n";
		$res .= date() . "\tmd5sum $LCY$folder1/$labelz.$ext.gz$N > $YW$folder1/.$labelz.$ext.gz.md5$N\n";
		$data->{$lbc}{lbcpreproc} .= "gzip -f $folder1/$labelz.$ext\n";
		$res .= date() . "\tmd5sum $LCY$folder1/$labelz.$ext.gz$N > $YW$folder1/.$labelz.$ext.gz.md5$N\n";
		$data->{$lbc}{lbcpreproc} .= "md5sum $folder1/$labelz.$ext.gz > $folder1/.$labelz.$ext.gz.md5\n";
		$ext .= ".gz";
	}
	$data->{$lbc}{lbcpreproc} =~ s/\/+/\//g;
#	$res .= "\t$data->{$lbc}{dedupe}\n";
#	$res .= "\t$data->{$lbc}{map}\n";
#	$res .= "\t$data->{$lbc}{peak}\n";
	$data->{$lbc}{dedupe} =~ s/DEDUPEINPUT/$resFile/g;
	$data->{$lbc}{dedupe} =~ s/DEDUPEOUTPUT/$resFile.rmdup.fq/g;
	$data->{$lbc}{dedupe} =~ s/DEDUPEOUTMD5/$folder1\/.$labelz.$ext.rmdup.fq.gz.md5/g;
	$data->{$lbc}{dedupe} =~ s/\/+/\//g;
	$res .= date() . "\t" . $data->{$lbc}{dedupe};
	$data->{$lbc}{map} =~ s/DEDUPEOUTPUT/$resFile.rmdup.fq/g;
	$data->{$lbc}{map} =~ s/FOOTLOOPOUTPUT/$resFile.rmdup.fq.gz/g;
	$data->{$lbc}{map} =~ s/\/+/\//g;
	$res .= date() . "\t" . $data->{$lbc}{map};
	$data->{$lbc}{peak} =~ s/FOOTLOOPOUTPUT/$resFile.rmdup.fq.gz/g;
	$data->{$lbc}{peak} =~ s/FOOTPEAKOUTPUT/$resFile.rmdup.fq.gz/g;
	$data->{$lbc}{peak} =~ s/\/+/\//g;
	$res .= date() . "\t" . $data->{$lbc}{peak};
	$res =~ s/\/+/\//g;
	$res .= "\n";
	$resFile =~ s/\/+/\//g;
	
	open (my $outSbatchDedupe, ">", "$resFile\_dedupe.sbatch") or DIELOG($outLog, "\n\nFailed to write sbatch $resFile.sbatch: $!\n");
	print $outSbatchDedupe "\#\!/bin/bash -l
\#SBATCH -p high --mem 48000

#Folder: $folder0
#File: $fileName0
#File2: $resFile
#uuid: $uuid
#date: $date
#run_script: $run_script

#preprocess
$data->{$lbc}{lbcpreproc}
#Dedupe
$data->{$lbc}{dedupe}

# map and peak
sbatch $resFile\_map.sbatch
";
	close $outSbatchDedupe;

	open (my $outSbatchMap, ">", "$resFile\_map.sbatch") or DIELOG($outLog, "\n\nFailed to write sbatch $resFile.sbatch: $!\n");
	print $outSbatchMap "\#\!/bin/bash -l
\#SBATCH -p high --mem 16000

#Folder: $folder0
#File: $fileName0
#File2: $resFile
#uuid: $uuid
#date: $date
#run_script: $run_script

#Map
$data->{$lbc}{map}
#Peak
$data->{$lbc}{peak}
";
	close $outSbatchMap;
	return ($res, $resFile);
}	


__END__
	my ($folder1, $fileName1) = getFilename($input1, "folder");
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		my @arr = split("\t", $line);
	}
	close $in1;
	
	
	open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
	close $out1;
	
	
	__END__
	0: 1
	1: PFC 8-1
	2: 1479
	3: 44
	4: 6.8
	5: BC1
	6: TCAGACGATGCGTCAT
	7: No transcription
	8: PFC 8
	9: CCCACTACGTGAACCATCACCC
	10: ATACGCAAACCGCCTCTCCCCG
