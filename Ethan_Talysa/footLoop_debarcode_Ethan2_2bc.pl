#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use File::Basename qw(dirname); use Cwd qw(abs_path);
use vars qw($opt_v $opt_b $opt_i $opt_d $opt_x $opt_2 $opt_w $opt_0 $opt_o);
getopts("vb:i:d:x:2:w:0o:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
}

use myFootLib; use FAlite; use footPeakAddon;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/../footLoop";
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
my ($barcodeFile2) = $opt_2 if defined $opt_2;
my $bc02;
my $usage = "-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N ${LGN}[-0: Dry run| -d: debug| -x fa.fai for bam]$N $CY-i <input1>$N $LPR-b <barcode.csv>$N

${LGN}Optional: 
-d <name for debug e.g. 203>
-x <fa.fai file if bam>$N\n";
#-2 <2nd barcode, 2 columns of ID and barcode sequence>

my $indexFile = $opt_x;

print "\n";
#print "$LRD ERROR$N -i not defined\n" if not defined $opt_i;
#print "$LRD ERROR$N -b not defined\n" if not defined $opt_b;

#die $usage unless defined $opt_i and -e $opt_i and defined $opt_b;
#die $usage if not defined $opt_i or not defined $opt_b;

die "$usage\n${LRD}ERROR$N: Does not exist -b ${LPR}opt_b$N!\n" if not defined $opt_b;
print "-b: $LPR$opt_b$N\n" if defined $opt_b;
die "$usage\n${LRD}ERROR$N: Does not exist -i $LCY$opt_i$N!\n" if defined $opt_i and not -e $opt_i;
print "-i: $LCY$opt_i$N\n" if defined $opt_i;
die "\n\n" if not defined $opt_0 and (not defined $opt_i or not defined $opt_b or (defined $opt_i and not -e $opt_i));
die "\n\n" if defined $opt_0 and (not defined $opt_b);

my $input1 = $opt_i;
my $debugname = $opt_d;
$debugname = "NA" if not defined $opt_d;
my $barcodeFile = $opt_b;
my @barcodeFiles = split(",", $barcodeFile);
my ($barcodeFileName) = getFilename($barcodeFile);
if (defined $opt_0 and not defined $opt_i) {
	system("touch dummy.temp") == 0 or die "Failed to create $LGN dummy.temp$N: $!\n";
	$opt_i = "dummy.temp";
	$input1 = $opt_i;
}

my ($inputName) = getFilename($input1);
if (defined $opt_o) {
	$inputName = $opt_o;
	$inputName =~ s/^\.+\/+//;
}

print "inputName=$YW$inputName$N\n";
die "Please supply with .fa.fai file since input1 is a bam/sam file!\n" if $input1 =~ /(.sam$|.bam$)/ and not defined $opt_x;
die "Please supply with .fa.fai file since input1 is a bam/sam file!\n" if $input1 =~ /(.sam$|.bam$)/ and defined $opt_x and not -e $opt_x;
print "-x: indexFile $opt_x\n" if defined $opt_x;
my $dbfolder = "$inputName\_barcode$barcodeFileName\_debarcode_result";
$dbfolder =~ s/\./_/g;
mkdir $dbfolder if not -d $dbfolder;
print "dbFolder=$YW$dbfolder$N\n";
my $outLogFile = "$inputName\_barcode$barcodeFileName\_logfile"; $outLogFile =~ s/\./_/g; $outLogFile = "$dbfolder/$outLogFile.txt";
my $outExtraFile = "$inputName\_barcode$barcodeFileName\_ExtraBC"; $outExtraFile =~ s/\./_/g; $outExtraFile = "$dbfolder/$outExtraFile.txt";
print "outlog File  : $LGN$outLogFile$N\n";
print "outExtra File: $LGN$outExtraFile$N\n";
open (my $outLog, ">", $outLogFile) or die "Cannot write to $LCY$outLogFile$N: $!\n";
open (my $outExtra, ">", $outExtraFile) or die "Cannot write to $LCY$outExtraFile$N: $!\n";
# parse barcode
my ($bc0, %out, %outcheck);
my @line = `cat @barcodeFiles`; my $linecount = -1;
print "\nOutput:\n";
my %bc;
foreach my $line (@line[0..@line-1]) {
	chomp($line);
	$linecount ++;
	next if $line =~ /^#/;
	next if $line =~ /^\n$/;
	next if $line !~ /[a-zA-Z0-9]+/;
	next if $line =~ /^(\#|No)/i or $line =~ /\ttarget/i;
	#my $delim = $line =~ /\t/ ? "\\t" : $line =~ /,/ ? "," : $line =~ /[ ]/ ? " " : die "Cannot determine delimiter at line=$LCY\n$line\n$N\n";
	#print "Delimiter = $delim\n";
	my @arr = split("[\t ]+", $line);
	my ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;
	my ($bcseq1, $bcseq2);
	if (@arr == 2) {
		($bc, $bcseq, $desc, $plasmid) = @arr;
	}
	else {
		($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#split(",", $line)  if $line !~ /\t/;
		if ($bcseq =~ /,/) {
			($bcseq1, $bcseq2) = split(",", $bcseq);
		}
#	   ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#osplit("\t", $line) if $line =~ /\t/;
		$desc =~ s/[ ]+//g;
		$plasmid =~ s/[ ]+//g;
	}
	$desc = "NA" if not defined $desc;
	$plasmid = "NA" if not defined $plasmid;
	$target = $bc if not defined $target;
	$desc =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	#$bcseq =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$bc =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$plasmid =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$target =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$bc = uc($bc);
	$plasmid = uc($plasmid);
	$desc = uc($desc);
	my $bcbefore = $bc;
	$bc = "bc$bc\_plasmid$plasmid\_desc$desc" if $plasmid ne "NA" or $desc ne "NA";#defined $plasmid and defined $desc;
	die "\n${LRD}FATAL ERROR$N: Multiple barcode name for the same barcode sequence.\n$LCY" . $bc0->{uc($bcseq)}{bc} . "$N\ncurrent=$LPR$bc$N\nbarcode=$LGN" . uc($bcseq) . "$N\n\n" if defined $bc0->{uc($bcseq)};
	if (@arr == 2) {
		$bc0->{uc($bcseq)}{plasmid} = $plasmid;
		$bc0->{uc($bcseq)}{target} = $target;
		$bc0->{uc($bcseq)}{total} = 0;
		$bc0->{uc($bcseq)}{line} = $linecount;
		$bc0->{uc($bcseq)}{bc} = $bc;
		#if (defined $bcseq2) {
		#	my $bc2;
		#	$bc0->{uc($bcseq2)}{plasmid} = $plasmid;
		#	$bc2->{uc($bcseq2)}{target} = $target;
		#	$bc2->{uc($bcseq2)}{total} = 0;
		#	$bc2->{uc($bcseq2)}{line} = $linecount;
		#	$bc2->{uc($bcseq2)}{bc} = $bc;
		#	$bc0->{uc($bcseq)}{bc2} = $bc2;
		#}
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc, bcseq=$LCY$bcseq$N\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N, bcseq=$LCY$bcseq$N\n";
	}
	else {
		$bc0->{uc($bcseq)}{plasmid} = $plasmid;
		$bc0->{uc($bcseq)}{target} = $target;
		$bc0->{uc($bcseq)}{bc} = $bc;
		$bc0->{uc($bcseq)}{total} = 0;
		$bc0->{uc($bcseq)}{line} = $linecount;
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc, bcseq=$LCY$bcseq$N\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N, bcseq=$LCY$bcseq$N\n";
		#if (defined $bcseq2) {
		#	my $bc2;
		#	$bc0->{uc($bcseq2)}{plasmid} = $plasmid;
		#	$bc2->{uc($bcseq2)}{target} = $target;
		#	$bc2->{uc($bcseq2)}{total} = 0;
		#	$bc2->{uc($bcseq2)}{line} = $linecount;
		#	$bc2->{uc($bcseq2)}{bc} = $bc;
		#	$bc0->{uc($bcseq)}{bc2} = $bc2;
		#	LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc, bcseq=$LCY$bcseq$N, bcseq2=$bcseq2\n","NA");
		#	print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N, bcseq=$LCY$bcseq$N\n";
		#}
		#else {
		#	LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc, bcseq=$LCY$bcseq$N\n","NA");
		#	print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N, bcseq=$LCY$bcseq$N\n";
		#}
	}
	#open ($out{$outBCFile}, ">", "$dbfolder/$inputName\_$bc.fastq") or die "Cannot write to $dbfolder/$inputName\_$bc.fastq: $!\n";
}
print $LGN . "\t" . (scalar(@line)) . "\t$dbfolder/UNKNOWN.fastq$N (no barcode)\n\n";

if (defined $opt_2) {
	my @line2 = `cat $barcodeFile2`; $linecount = -1;
	foreach my $line2 (@line2[0..@line2-1]) {
		chomp($line2);
		$linecount ++;
		next if $line2 =~ /^\n$/;
		next if $line2 !~ /[a-zA-Z0-9]+/;
		my ($target2, $bcseq2, $plasmid2) = split("\t", $line2);
		$bc02->{uc($bcseq2)}{plasmid} = $plasmid2;
		$bc02->{uc($bcseq2)}{target} = $target2;
		$bc02->{uc($bcseq2)}{bc} = uc("$target2");
		$bc02->{uc($bcseq2)}{total} = 0;
		$bc02->{uc($bcseq2)}{line} = $linecount;
		foreach my $bcseq (sort keys %{$bc0}) {
			my $plasmid = $bc0->{$bcseq}{plasmid};
			next if $plasmid ne $plasmid2;
			if (defined $bc0->{$bcseq}{bc2}) {
				if (defined $bc0->{$bcseq}{bc2}{uc($bcseq2)}) {
					die "Alreadyd defined bcseq2=$LCY$bcseq2$N for bc=$bc0->{$bcseq}{bc} (bcseq1=$bcseq)\n
target=$bc0->{$bcseq}{bc2}{uc($bcseq2)}{target}
curr target=$target2";
				}
			}
			$bc0->{$bcseq}{bc2}{uc($bcseq2)}{plasmid} = $plasmid2;
			$bc0->{$bcseq}{bc2}{uc($bcseq2)}{target} = $target2;
			$bc0->{$bcseq}{bc2}{uc($bcseq2)}{bc} = uc($target2);
			$bc0->{$bcseq}{bc2}{uc($bcseq2)}{total} = 0;
			$bc0->{$bcseq}{bc2}{uc($bcseq2)}{line} = $linecount;
			#print "bc1=$bc0->{$bcseq}{bc}, bcseq1=$bcseq, plasmid1=$bc0->{$bcseq}{plasmid}\nbc2=$bc02->{uc($bcseq2)}{bc}, bcseq2=$bcseq2, plasmid2=$bc02->{uc($bcseq2)}{plasmid}\n\n";
		}
	}
}



if (defined $opt_0) {
	print "\nDry run ${LGN}SUCCESS$N\n\n";
	exit;
}

open (my $outnobc, ">", "$dbfolder/$inputName\_UNKNOWN.fastq") or die "Cannot write to $dbfolder/$inputName\_UNKNOWN.fastq: $!\n";
my ($seq_nobc, $seq_dblbc, $seq_totbc, $readCount) = (0,0,0,0);
my ($data, $bad, @len, $in1, %count);
my ($folder1, $fileName1) = getFilename($input1, "folder");
print "${LGN}1. -b $LCY$barcodeFile$N was parsed successfully!\n$N\n";
# parse fastq
my $isDir = -d $input1 ? "(is a directory!)" : -e $input1 ? "" : "($input1 does not exist!)";
die "$LRD FATAL ERROR$N: -i $LCY$input1$N is not a fastq file! $LRD$isDir$N\n" if -d $input1 or not -e $input1;
if ($input1 =~ /(.fq.gz$|.fq$|.fastq$|.fastq.gz$)/) {
	open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.gz$/i;
	open ($in1, "zcat < $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.gz$/i;
	while (my $line = <$in1>) { chomp($line);
		# 1. name: m171023_184836_42145_c101402632550000001823289101251862_s1_p0/0/ccs
		$readCount ++;
		my ($readName) = $line;
		$line = <$in1>; chomp($line);
		# 2. seqs: ACTGTCG....
		my ($seq) = $line;
		my $len = int(length($seq));
		
		$line = <$in1>; chomp($line);
		# 3. junk (+)

		$line = <$in1>; chomp($line);
		# 4. Quality
		my $qua = $line;

		my ($bc, $bcs, $seq2, $qua2, $bc3, $goodtype);
		next if defined $opt_d and $readName !~ /$debugname/;
		($bc, $bcs, $seq2, $qua2, $bc3, $goodtype) = infer_barcode($bc0, $seq, $qua, $readName, $data, $bc02);
	#return ($barcode, $barcode2, $seq, $qua, $barcode3, $goodtype);
		$bc3 = "UNKNOWN" if not defined $bc3;
		$bc = "UNKNOWN" if not defined $bc;
		print $outLog "$readCount: $LGN$readName$N, bc=$bc, bc3=$bc3\n";
		print $outnobc "$readName\n$seq\n\+\n$qua\n" if not defined $seq2;
		print $outLog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
		print $outExtra "$readName\t$bc\t$bc3\t$goodtype\n";
		my $outBC = "$inputName\_$bc"; $outBC =~ s/\./_/g; 
		my $outBCFile = "$dbfolder/$outBC.fq.gz";

		$count{$outBC}{total} ++;

		if (not defined $outcheck{$outBCFile}) {
			open ($out{$outBCFile}, "| gzip > $outBCFile") or die "Cannot write to $LCY$outBCFile$N: $!\n";
			$outcheck{$outBCFile} = 1;
		}
		print {$out{$outBCFile}} "$readName\n$seq2\n\+\n$qua2\n";# if defined $seq2;

		print date() . "\t$YW$input1$N: Done $LGN$readCount$N\n" if $readCount % 500 == 0;
		last if $readCount > 6;
	}
	close $in1;
}
elsif ($input1 =~ /(.sam$|.bam$)/) {
	open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.bam$/i;
	open ($in1, "samtools view $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.bam$/i;
	while (my $line = <$in1>) { chomp($line);
		my @arr = split("\t", $line);
		next if @arr < 10;
		$readCount ++;
		# 1. name: m171023_184836_42145_c101402632550000001823289101251862_s1_p0/0/ccs
		# 2. seqs: ACTGTCG....
		my ($readName, $plasmid, $seq, $qua) = ($arr[0], $arr[2], $arr[9], $arr[10]);
		my $len = int(length($seq));

		my ($bc, $bcs, $seq2, $qua2, $bc3, $goodtype);
		next if defined $opt_d and $readName !~ /$debugname/; #DEBUG
		($bc, $bcs, $seq2, $qua2, $bc3, $goodtype) = infer_barcode($bc0, $seq, $qua, $readName, $data, $bc02);
#	return ($barcode, $barcode2, $seq, $qua, $barcode3, $goodtype);
		#return ($barcode, $barcode2, $seq, $qua, $barcode3);
		$bc3 = "UNKNOWN" if not defined $bc3;
		$bc  = "UNKNOWN" if not defined $bc;
		
		#if (not defined $out{$outBCFile}) {
		my $outBC = "$inputName\_$bc"; $outBC =~ s/\./_/g; 
		$outBC =~ s/plasmid.+_desc/plasmid$plasmid\_desc/ if $outBC =~ /plasmid.+_desc/;
		$outBC =~ s/_UNKNOWN/_plasmid$plasmid\_descUNKNOWN/ if $outBC =~ /UNKNOWN/;
		$outBC =~ s/_bam_/_/;


		$count{$outBC}{total} ++;

		print $outLog "$readCount: $LGN$readName$N, BC=$outBC, bc=$bc\n";
		print $outnobc "$readName\n$seq\n\+\n$qua\n" if not defined $seq2;
		print $outLog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
		print $outExtra "$readName\t$outBC\t$bc\t$goodtype\n";

		my $outBCFile = "$dbfolder/$outBC.fq";
		if (not defined $outcheck{$outBCFile}) {
			$outcheck{$outBCFile} = 1;
			LOG($outLog, "OUTBCFILE=$LCY$outBCFile$N\n");
			open ($out{$outBCFile}, "> $outBCFile") or die "Cannot write to $LCY$outBCFile$N: $!\n";
			#open ($out{$outBCFile}, "| samtools view -bS -t $indexFile - > $outBCFile") or die "Cannot write to $LCY$outBCFile$N: $!\n";
		}
		print {$out{$outBCFile}} "\@$readName\n$seq\n+\n$qua\n";
		#print {$out{$outBCFile}} "$line\n";

		print date() . "\t$YW$input1$N: Done $LGN$readCount$N\n" if $readCount % 500 == 0;
		#last if $readCount > 6;
	}
	close $in1;
}

foreach my $outBCFile (sort keys %out) {
	close ($out{$outBCFile});
}

close $outExtra;
$seq_totbc = $readCount;
print "${LGN}2. -i $LCY$input1$N was parsed successfully!\n$N\n";
my $double_bc = (keys %{$bad});
my $seq_goodbc = $seq_totbc - $seq_nobc - $double_bc;
my $seq_goodbcPerc = $seq_totbc == 0 ? 0 : int($seq_goodbc / $seq_totbc * 1000)/10;
print STDERR "

Total seq : $seq_totbc
Has bc    : $seq_goodbc ($seq_goodbcPerc \%)
Double bc : $double_bc
No bc     : $seq_nobc

";
print $outLog "

Total seq : $seq_totbc
Has bc    : $seq_goodbc ($seq_goodbcPerc \%)
Double bc : $double_bc
No bc     : $seq_nobc

";

foreach my $outBC (sort {$count{$b} <=> $count{$a}} keys %count) {
	print $outLog "\t$outBC\t$count{$outBC}{total}\n";
	print STDERR "\t$outBC\t$count{$outBC}{total}\n";
}
#foreach my $bcseq (sort {$bc0->{$a}{line} <=> $bc0->{$b}{line}} keys %{$bc0}) {
#	print $outLog "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
#	print STDERR "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
#}
sub myrevcomp {
	my ($seq) = @_;
	$seq =~ tr/ACGT/TGCA/;
	$seq = rev2($seq);
	return($seq);
}

sub rev2 {
	my ($seq) = @_;
	my $seq0 = $seq;
	if ($seq =~ /\.\{\d+,?\d*\}/) {
		my ($seq1, $re, $seq2) = $seq =~ /^(.*)(\.\{\d+,?\d*\})(.*)$/;
		$seq1 = "" if not defined $seq1;
		$seq2 = "" if not defined $seq2;
		$seq1 = reverse($seq1);
		$seq2 = reverse($seq2);
		return($seq2 . $re . $seq1);
	}
	elsif ($seq =~ /^\^(.+)\.\+(.+)\$$/) {
		my ($seq1, $re, $seq2) = $seq =~ /^\^(.*)(\.\+)(.*)\$$/;
		$seq1 = "" if not defined $seq1;
		$seq2 = "" if not defined $seq2;
		$seq1 = reverse($seq1);
		$seq2 = reverse($seq2);
		return("^$seq2$re$seq1\$");
	}
	else {
		die "seq0=\n$seq0\n$seq\n\n" if $seq !~ /^[NACTG]+$/;
		return(reverse($seq));
	}
}

sub testbc {
	my ($bc, $bcseq, $nameshort, $seq) = @_;
	my $bcseqrev = $bcseq; 
	$bcseqrev =~ tr/ACTG/TGAC/;
	$bcseqrev = rev2($bcseqrev);
				
	print $outLog "\n${YW}1. TEST FORWARD $bcseq$N\n";
	my $muscleinput = ">$bc\n$bcseq\n>$nameshort\n$seq";
	my @res = muscle($muscleinput);
	#print "RES:\n$LGN" . join("", @res) . "\n$N\n";
	#print "muscleinput:$LPR$muscleinput$N\n";
	my ($matchperc, $pos) = parse_muscle(\@res, $bcseq, $seq);
	my $rez2 = "";
	for (my $r = 0; $r < @res; $r++) {
		chomp($res[$r]);
		if ($res[$r] =~ /^>/) {
			$res[$r] .= "MUHSPACE";
		}
		$rez2 .= $res[$r];
	}
	$rez2 =~ s/>/\n/g;
	$rez2 =~ s/^\n$//g;
	$rez2 =~ s/MUHSPACE/\n/g;
	return($matchperc, $pos, $rez2);
}

sub testbc2 {
	my ($testbc, $testbcseq, $nameshort, $seq) = @_;
				
	#print $outLog "\n${YW}1. TEST FORWARD $testbcseq$N\n";
	#print "\n${YW}1. TEST FORWARD $testbcseq$N\n";
	my $muscleinput = ">$testbc\n$testbcseq\n>$nameshort\n$seq";
	my @res = muscle($muscleinput);
	#print $outLog "RES:\n$LGN" . join("", @res) . "\n$N\n";
	#print $outLog "muscleinput:$LPR$muscleinput$N\n";
	#print "RES:\n$LGN" . join("", @res) . "\n$N\n";
	#print "muscleinput:$LPR$muscleinput$N\n";
	my ($matchperc, $pos) = parse_muscle(\@res, $testbcseq, $seq);
	my $rez2 = "";
	for (my $r = 0; $r < @res; $r++) {
		chomp($res[$r]);
		if ($res[$r] =~ /^>/) {
			$res[$r] .= "MUHSPACE";
		}
		$rez2 .= $res[$r];
	}
	$rez2 =~ s/>/\n/g;
	$rez2 =~ s/^\n$//g;
	$rez2 =~ s/MUHSPACE/\n/g;
	return($matchperc, $pos, $rez2, $testbcseq);
}

sub parse_bcseq {
	my ($bcseqs) = @_;
	my %bc;
	my @bcseqs = split(",", $bcseqs);
	for (my $i = 0; $i < @bcseqs; $i++) {
		my $bcseq = $bcseqs[$i];
		if ($bcseq =~ /^(BEG|END)/) {
			my ($bctype, $number, $seq) = $bcseq =~ /^(BEG|MID|END)(\d+)=(.+)$/;
			$bctype = $bctype eq "BEG" ? "1BEG" : $bctype eq "MID" ? "2MID" : $bctype eq "END" ? "3END" : die "bctype ($bctype) has to be BEG/MID/END!\n";
			$bc{sep}{$bctype}{$number} = $seq;
		}
	}
	foreach my $bctype (sort keys %{$bc{sep}}) {
		foreach my $number (sort {$a <=> $b} keys %{$bc{sep}{$bctype}}) {
			$bc{all}{$bctype} .= $bc{sep}{$bctype}{$number};
		}
	}
	#foreach my $bctype (sort keys %{$bc{all}}) {
	#	print "$bctype = $bc{all}{$bctype}\n";
	#}
	return \%bc;
}



sub infer_barcode {
	my ($bc1, $seq, $qua, $name, $data, $bc02) = @_;
	LOG($outLog, "\n---------------------------------------\n${YW}$name$N\n","NA");
	my $barcode; my $barcode2; my $barcode3; my $type;
	my $seq1; my $seq2;
	my $seq1chunk; my $seq2chunk;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	my %good;
	$good{highest} = 0;
	my $pos = -1;
	my $matchperc = 0;
	my $goodtype = "BAD";
	$barcode3 = "$LRD$matchperc$N\t$LPR$pos$N";

	foreach my $bcseqorig (sort {$bc1->{$a}{bc} cmp $bc1->{$b}{bc}} keys %{$bc1}) {
		my $bcseq = $bcseqorig;
		my $bc = $bc1->{$bcseq}{bc};
		my $nameshort = $name; $nameshort =~ s/^\@//;

		if ($bcseq =~ /[\=,]/) {
			my $bad = 0;
			my $bclong = parse_bcseq($bcseq);
			die "bc=$bc: at infer_bc, undef bclong from bcseq=$LGN$bcseq$N\n" if not defined $bclong;
			my %bclong = %{$bclong};
			my $bcfwall = "";
			my $bcrvall = "";
			my $bcstrand = "FW";
			foreach my $bctype (sort keys %{$bclong{all}}) {
				my $bcfw = $bclong{all}{$bctype};
				my $bcrv = myrevcomp($bcfw);
				if ($bctype eq "1BEG") {
					#print " 1BEG: test\n$LGN^$bcfw$N\n$LCY$bcrv\$$N\n";
					$bad ++ if $seq !~ /(^$bcfw|$bcrv$)/;
					my $bcstrand2 = $seq =~ /($bcrv$)/ ? "RV" : "FW";
					$bcstrand = $bcstrand2 if not defined $bcstrand;
					$bad ++ if defined $bcstrand and $bcstrand ne $bcstrand2;
				}
				elsif ($bctype eq "2MID") {
					$bad ++ if $seq !~ /($bcfw|$bcrv)/;
					my $bcstrand2 = $seq =~ /($bcrv)/ ? "RV" : "FW";
					$bcstrand = $bcstrand2 if not defined $bcstrand;
					$bad ++ if defined $bcstrand and $bcstrand ne $bcstrand2;
				}
				elsif ($bctype eq "3END") {
					#print " 3END: test\n$LGN^$bcfw$N\n$LCY$bcrv\$$N\n";
					$bad ++ if $seq !~ /($bcfw$|^$bcrv)/;
					my $bcstrand2 = $seq =~ /(^$bcrv)/ ? "RV" : "FW";
					$bcstrand = $bcstrand2 if not defined $bcstrand;
					$bad ++ if defined $bcstrand and $bcstrand ne $bcstrand2;
				}
				$bcfwall .= $bcfwall eq "" ? "$bcfw" : ",$bcfw";
				$bcrvall .= $bcrvall eq "" ? "$bcrv" : ",$bcrv";
			}
			if ($bad == 0) {
				$pos = 0;
				$matchperc = 100;
				$goodtype = "${LGN}PERFECT_FW_$bcfwall$N" if $matchperc == 100 and $bcstrand eq "FW";
				$goodtype = "${LGN}PERFECT_RV_$bcrvall$N" if $matchperc == 100 and $bcstrand eq "RV";
				$barcode3 = "$LGN$matchperc$N\t$LPR$pos$N";
				$good{highest} = $matchperc;
				undef $good{bc};
				$good{bc}{$bc} = $bcseq;
				#print "$goodtype\n";
				last;
			}
			else {
				my @split = split("[,\=]", $bcseq);
				
				#print "bad bc $LCY$bcseq$N\n";
				#print join("\n", @split) . "\n\n";
			}
		}
		else {
			my $bcseqrev = $bcseq; 
			$bcseqrev =~ tr/ACTG/TGAC/;
			$bcseqrev = rev2($bcseqrev);
			
			if ($seq =~ /($bcseq|$bcseqrev)/) {
				($pos) = $seq =~ /^(.*)($bcseq|$bcseqrev)/;
				$pos = defined $pos ? length($pos) : -1;
				$matchperc = 100;
				$goodtype = "${LGN}PERFECT_FW_$bcseq$N" if $matchperc == 100 and $seq =~ /$bcseq/;
				$goodtype = "${LGN}PERFECT_RV_$bcseqrev$N" if $matchperc == 100 and $seq =~ /$bcseqrev/;
				$barcode3 = "$LGN$matchperc$N\t$LPR$pos$N";
				$good{highest} = $matchperc;
				undef $good{bc};
				$good{bc}{$bc} = $bcseq;
				last;
			}
		}
	}

	#print "matchperc=$matchperc\n";
	if ($matchperc < 100) {
		$goodtype = "BAD";

		foreach my $bcseqorig (sort {$bc1->{$a}{bc} cmp $bc1->{$b}{bc}} keys %{$bc1}) {

			my $bc = $bc1->{$bcseqorig}{bc};
			my $nameshort = $name; $nameshort =~ s/^\@//;
			my $muscleinput;
			my @res;
			my ($matchperc, $pos, $rez2, $bcseq, $bcseqrev, $beststrand);

			if ($bcseqorig =~ /[\=,]/) {
				my $bad = 0;
				my $bclong = parse_bcseq($bcseqorig);
				die "bc=$bc: at infer_bc, undef bclong from bcseq=$LGN$bcseqorig$N\n" if not defined $bclong;
				my %bclong = %{$bclong};
				my $bcfwall = "";
				my $bcrvall = "";
				my $bcstrand = "FW";
				my $matchpercfw0 = 0;
				my $matchpercrv0 = 0;
				my $totalmatch = 0;
				my $posfw0 = -1;
				my $posrv0 = -1;
				my $rez2fw0 = "";
				my $rez2rv0 = "";
				my $bcseqfw0 = "";
				my $bcseqrv0 = "";
				#test FW
				foreach my $bctype (sort keys %{$bclong{all}}) {
					my $bcfw = $bclong{all}{$bctype};
					my $bcrv = myrevcomp($bcfw);
					my ($matchpercfw, $posfw, $rez2fw, $bcseqfw) = testbc2($bc, $bcfw, $nameshort, $seq);
					my ($matchpercrv, $posrv, $rez2rv, $bcseqrv) = testbc2($bc, $bcrv, $nameshort, $seq);
					LOG($outLog, "- bc=$YW$bc$N, $bctype:$LPR FW$N: $LGN$matchpercfw$N % match at $LGN$posfw$N/$LCY$seq_length$N,$LPR RV$N: $LGN$matchpercrv$N % match at $LGN$posrv$N/$LCY$seq_length$N\n","NA");
					$posfw0 = $posfw if ($posfw0 eq -1 or ($posfw0 ne -1 and $posfw0 < $posfw));
					$posrv0 = $posrv if ($posrv0 eq -1 or ($posrv0 ne -1 and $posrv0 < $posrv));
					$matchpercfw0 += $matchpercfw;
					$matchpercrv0 += $matchpercrv;
					$rez2fw0 .= $rez2fw0 eq "" ? $rez2fw : "\n$rez2fw";
					$rez2rv0 .= $rez2rv0 eq "" ? $rez2rv : "\n$rez2rv";
					$bcseqfw0 .= $bcseqfw0 eq "" ? $bcseqfw : ".+$bcseqfw";
					$bcseqrv0 .= $bcseqrv0 eq "" ? $bcseqrv : ".+$bcseqrv";
					$totalmatch ++;
					$bcfwall .= $bcfwall eq "" ? "$bcfw" : ",$bcfw";
					$bcrvall .= $bcrvall eq "" ? "$bcrv" : ",$bcrv";
				}
				
				$matchpercfw0 = $totalmatch eq 0 ? 0 : int(10*$matchpercfw0 / $totalmatch+0.5)/10;
				$matchpercrv0 = $totalmatch eq 0 ? 0 : int(10*$matchpercrv0 / $totalmatch+0.5)/10;
				if ($matchpercrv0 <= $matchpercfw0) {
					$matchperc = $matchpercfw0;
					$pos = $posfw0;
					$rez2 = $rez2fw0;
					$bcseq = $bcseqfw0;
					$bcseqrev = $bcseqrv0;
					$beststrand = "FW";
				}
				else {
					$matchperc = $matchpercrv0;
					$pos = $posrv0;
					$rez2 = $rez2rv0;
					$bcseq = $bcseqfw0;
					$bcseqrev = $bcseqrv0;
					$beststrand = "RV";
				}
				LOG($outLog, "$name: matchperc=$YW$matchperc$N pos=$LPR$pos$N bcseq=$LGN$bcseq $LCY$bcseqrev$N\n\n","NA");
			}
			else {
				my ($matchpercfw, $posfw, $rez2fw) = testbc($bc, $bcseqorig, $nameshort, $seq);
				my ($matchpercrv, $posrv, $rez2rv) = testbc($bc, myrevcomp($bcseqorig), $nameshort, $seq);
				($matchperc, $pos, $rez2) = $matchpercfw >= $matchpercrv ? ($matchperc, $posfw, $rez2fw) : ($matchpercrv, $posrv, $rez2rv);
				($bcseq, $bcseqrev) = ($bcseqorig, myrevcomp($bcseqorig));
			}
				
			LOG($outLog, "$rez2\n","NA");
			#print join("", @res) . "\n\n";
			if ($matchperc > $good{highest}) {
				$barcode3 = "$LGN$matchperc$N\t$LPR$pos$N";
				$barcode3 = "$LRD$matchperc$N\t$LPR$pos$N" if $matchperc < 50;
				$barcode3 = "$YW$matchperc$N\t$LPR$pos$N" if $matchperc >= 50 and $matchperc < 86;
				$good{highest} = $matchperc;
				undef $good{bc};
				$goodtype = $beststrand eq "FW" ? "${LRD}FW_$bcseq$N" : "${LRD}RV_$bcseqrev$N";
				$good{bc}{$bc} = $bcseq;
				last if $matchperc >= 99;
			}
		}
	}
	
	#$good2{highest} = 0;
	foreach my $bc (sort keys %{$good{bc}}) {
		my $bcseq = $good{bc}{$bc};
		#print $outLog "  BC=$YW$bc$N\t${LGN}GOOD$N=$good{bc}{$bc}$N highest=$LGN$good{highest}$N\n" if (keys %{$good{bc}}) == 1;
		#print $outLog "  BC=$YW$bc$N\t${LCY}DOUBLE BC$N=$good{bc}{$bc}$N highest=$LGN$good{highest}$N\n" if (keys %{$good{bc}}) > 1;
		LOG($outLog, "\n\n  BC=$YW$bc$N\t${LGN}GOOD$N=$good{bc}{$bc}$N ($barcode3,$goodtype) highest=$LGN$good{highest}$N\n") if (keys %{$good{bc}}) == 1;
		LOG($outLog, "\n\n  BC=$YW$bc$N\t${LCY}DOUBLE BC$N=$good{bc}{$bc}$N ($barcode3,$goodtype) highest=$LGN$good{highest}$N\n") if (keys %{$good{bc}}) > 1;
		#print "\n\n  BC=$YW$bc$N\t${LGN}GOOD$N=$good{bc}{$bc}$N ($barcode3,$goodtype) highest=$LGN$good{highest}$N\n" if (keys %{$good{bc}}) == 1;
		#print "\n\n  BC=$YW$bc$N\t${LCY}DOUBLE BC$N=$good{bc}{$bc}$N ($barcode3,$goodtype) highest=$LGN$good{highest}$N\n" if (keys %{$good{bc}}) > 1;
	}
	if ($good{highest} == 100 and (keys %{$good{bc}}) == 1) {
		foreach my $bc (keys %{$good{bc}}) {
			my $bcseq = $good{bc}{$bc};
			$barcode = $bc;
			$barcode2 = $bc;
			$bc0->{$bcseq}{total} ++;
		}
	}
	elsif ($good{highest} >= 86) {
		if ((keys %{$good{bc}}) == 1) {
			foreach my $bc (keys %{$good{bc}}) {
				my $bcseq = $good{bc}{$bc};
				$barcode = $bc;
				$barcode2 = $bc;
				$bc0->{$bcseq}{total} ++;
			}
		}
		else {
			$bad->{$name} ++;
			undef $barcode;
			undef $barcode2;
			undef $barcode3;
		}
	}
#foreach my $bcseq (sort {$bc0->{$a}{line} <=> $bc0->{$b}{line}} keys %{$bc0}) {
#	print $outLog "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
#	print STDERR "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
#}
	elsif ($good{highest} < 86) {
		#undef $barcode;
		undef $barcode2;
		#undef $barcode3;
		$seq_nobc ++;
	}
#	return ($barcode, $barcode2, $seq, $qua, $barcode3);
	#print "\n\nbarcode=$barcode\nbarcode2=$barcode2\nbarcode3=$barcode3\ngoodtype=$goodtype\n";
	return ($barcode, $barcode2, $seq, $qua, $barcode3, $goodtype);
#	($bc, $bcs, $seq2, $qua2, $data, $bc3) = infer_barcode($bc0, $seq, $qua, $readName, $data, $bc02);
}

sub parse_muscle {
	my ($res, $inputseq1, $inputseq2) = @_;
	my @res = @{$res};
	my ($name, $seq);
	my @name; my @seq;
	for (my $i = 0; $i < @res; $i++) {
#DEBUGDIE 
#die "res=$res[$i]\n";
		my $line = $res[$i];
		chomp($line);
		if ($line =~ /^>/) {
			if (defined $name) {
				push(@name, $name);
				push(@seq, $seq);
			}
			$name = $line;
			$seq = "";
		}
		else {
			$seq .= $line;
		}
	}
	if (defined $name) {
		push(@name, $name);
		push(@seq, $seq);
	}
	my $length = length($inputseq1) < length($inputseq2) ? length($inputseq1) : length($inputseq2);
	my $name1 = $name[0];
	my $name2 = $name[1];
	my $seq1 = $seq[0];
	my $seq2 = $seq[1];
#	die "
#inputseq1=$LCY$inputseq1$N
#inputseq2=$LCY$inputseq2$N
#name1=$LCY$name1$N
#name2=$LCY$name2$N
#seq1=$LCY$seq1$N
#seq2=$LCY$seq2$N
#";
	my ($seq1beg) = $seq1 =~ /^([\-]+)[A-Za-z0-9]/; $seq1beg = "" if not defined $seq1beg;
	my ($seq2beg) = $seq2 =~ /^([\-]+)[A-Za-z0-9]/; $seq2beg = "" if not defined $seq2beg;
	my ($seq1end) = $seq1 =~ /[A-Za-z0-9]([\-]+)$/; $seq1end = "" if not defined $seq1end;
	my ($seq2end) = $seq2 =~ /[A-Za-z0-9]([\-]+)$/; $seq2end = "" if not defined $seq2end;
	my $seq1beglen = length($seq1beg);
	my $seq2beglen = length($seq2beg);
	my $seq1endlen = length($seq1end);
	my $seq2endlen = length($seq2end);
	my $pos = 0;
#PRINT HERE
	#print "\n\n$LGN$name1$N\n$seq1\n$LCY$name2$N\n$seq2\n";
#PRINT HERE
	if (length($seq1beg) > 0) {
		my ($seq2a, $seq2b) = $seq2 =~ /^(.{$seq1beglen})(.+)$/;
		($pos) = $seq2a =~ tr/ACTGN/ACTGN/;
		$seq1 =~ s/^.{$seq1beglen}(.+)$/$1/;
		$seq2 =~ s/^.{$seq1beglen}(.+)$/$1/;	
#PRINT HERE
	#	print "\n${LCY}!POS=$N$LCY$pos$N\n$seq2a\n$seq2b\n";
#PRINT HERE
	}
	if (length($seq1end) > 0) {
		$seq1 =~ s/^(.+).{$seq1endlen}$/$1/;
		$seq2 =~ s/^(.+).{$seq1endlen}$/$1/;
	}
#	if (length($seq1beg) > length($seq2beg)) {
#		$seq1 =~ s/^.{$seq1beglen}(.+)$/$1/;
#		$seq2 =~ s/^.{$seq1beglen}(.+)$/$1/;
#	}
#	else {
#		$seq1 =~ s/^.{$seq2beglen}(.+)$/$1/;
#		$seq2 =~ s/^.{$seq2beglen}(.+)$/$1/;
#	}

#	if (length($seq1end) > length($seq2end)) {
#		$seq1 =~ s/^(.+).{$seq1endlen}$/$1/;
#		$seq2 =~ s/^(.+).{$seq1endlen}$/$1/;
#	}
#	else {
#		$seq1 =~ s/^(.+).{$seq2endlen}$/$1/;
#		$seq2 =~ s/^(.+).{$seq2endlen}$/$1/;
#	}
	#print "\n\n$LGN$name1$N\n$seq1\n$LCY$name2$N\n$seq2\n";
	
	my @seq1 = split("", $seq1);
	my @seq2 = split("", $seq2);
	my $beg = 0;
	my $match = 0; my $total = 0;
	for (my $i = 0; $i < @seq1; $i++) {
		$match ++ if $seq1[$i] eq $seq2[$i] and $seq1[$i] ne "-" and $seq1[$i] ne "N";
		$total ++ if $seq1[$i] ne "N";
	}
	my $matchperc = $total == 0 ? 0 : int($match / $total * 1000+0.5)/10;
#PRINT HERE
	#print "match = $match/$total ($matchperc %%)\n";
	#print "\n\n$LGN$name1$N\n$seq1\n$LCY$name2$N\n$seq2\n";
#PRINT HERE
	#die if length($seq1beg) > 0;
	return ($matchperc, $pos);
}
#		print "$i $LGN$res[$i]$N\n";#$name1,$seq1\n$name2,$seq2\n$name3,$seq3\n";
#	foreach my $name (sort {$res2{$a}{order} <=> $res2{$b}{order}} keys %res2) {
#		
#		print "$name,$res2{$name}{seq}\n";
#	}
#	die;
#}
#
#		chomp($res[$i]); $res[$i] =~ s/[ \t]+/\t/;
#		my ($name1, $seq1) = $res[$i] =~ /^(.+)\t(.+)$/;
#
#		$i ++; chomp($res[$i]); $res[$i] =~ s/[ \t]+/\t/;
#		my ($name2, $seq2) = $res[$i] =~ /^(.+)\t(.+)$/;
#
#		$i ++; chomp($res[$i]); $res[$i] =~ s/[ \t]+/\t/;
#		my ($name3, $seq3) = $res[$i] =~ /^(.+)\t(.+)$/;
#
#		die;
#	}
#
#}
sub muscle {
   my ($cmd) = @_;
   return(`echo '$cmd' | clustalo -i - 2> /dev/null`);
}

=comment
	for (my $l = 0; $l < 4; $l++) {
		foreach my $bcseq (sort {$bc1->{$a}{bc} cmp $bc1->{$b}{bc}} keys %{$bc1}) {
			my $bcseqchunk = substr($bcseq, $l, length($bcseq)-$l+1);
			my $length = length($bcseqchunk);
			my $bcseqchunkrev = $bcseqchunk;
			$bcseqchunkrev =~ tr/ACGT/TGCA/;
			$bcseqchunkrev = rev2($bcseqchunkrev);
			#if (($seq =~ /^(.){0,10}($bcseqchunk|$bcseqchunkrev)/ or $seq =~ /($bcseqchunk|$bcseqchunkrev).{0,10}$/)) {
			if ($seq =~ /($bcseqchunk|$bcseqchunkrev)/) {# or $seq =~ /($bcseqchunk|$bcseqchunkrev)/)) {
				if (defined $bc1->{$bcseq}{bc2}) {
					($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) = get_second_barcode($seq, $bc1, $bcseq, $name, $bad, $data, $bc02);
				}
				if (not defined $barcode) {
					my $pos = 0 if $seq =~ /^($bcseqchunk|$bcseqchunkrev)/;
					   $pos = $seq_length if $seq =~ /($bcseqchunk|$bcseqchunkrev)$/;
					my ($pos1, $pos2) = $seq =~ /^(.+)($bcseqchunk|$bcseqchunkrev)(.+)$/ if not defined $pos and $seq =~ /($bcseqchunk|$bcseqchunkrev)/;
						$pos = defined $pos ? $pos : length($pos1) <= 10 ? length($pos1) : length($pos2);
					my $type = $pos < 0.5 * $seq_length ? 1 : 2;
					die "$name: barcode=$barcode, bcseq=$bcseq ($bcseqchunk), pos1=$pos1 pos2=$pos2\n" if not defined $pos;
					$barcode2 = $seq =~ /$bcseqchunk/ ? "$type,$bcseqchunk" : "$type,$bcseqchunkrev";
					if (defined $barcode) {
						$bad->{$name} = 1;
						$seq2 = $bcseq;
						$seq2chunk = $bcseqchunk;
						print $outLog "$LCY$name\tDBL\t$pos/$seq_length\tprev=$barcode\tcurr=$bc1->{$bcseq}{bc}\tprev=$seq1\tcurr=$seq2$N\n";
					}
					else {
						$barcode = $bc1->{$bcseq}{bc};
						$seq1 = $bcseq;
						$seq1chunk = $bcseqchunk;
						$data->{$name} = $barcode;
						$bc1->{$bcseq}{total} ++;
						print $outLog "$LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n";
					}
				}
			}
			die if $l > 3;
		}
		$seq_totbc ++ if $l == 0;
		if (defined $barcode and not defined $bad->{$name}) {
			$barcode3 = $barcode if not defined $barcode3;
			my $bcseq = $seq1;
			my $bcseqchunk    = $seq1chunk;
			my $bcseqchunkrev = $seq1chunk;
			$bcseqchunkrev =~ tr/ACGT/TGCA/;
			$bcseqchunkrev = rev2($bcseqchunkrev);
			#my ($bad1, $good1) = $seq =~ /^(.{0,14}$bcseqchunk|.{0,14}$bcseqchunkrev)(.+)$/ if $seq =~ /^(.{0,14}$bcseqchunk|.{0,14}$bcseqchunkrev)/;
			#my ($good2, $bad2) = $seq =~ /^(.+)($bcseqchunk.{0,14}|$bcseqchunkrev.{0,14})$/ if $seq =~ /^(.+)($bcseqchunk.{0,14}|$bcseqchunkrev.{0,14})$/;
			#my ($bad1, $good1, $good2, $bad2);
			#if ($seq =~ /($bcseqchunk|$bcseqchunkrev)/) {
			#	($bad1, $good1) = $seq =~ /($bcseqchunk|$bcseqchunkrev)/;
			#}
			#else {
			#	($good2, $bad2) = $seq =~ /($bcseqchunk|$bcseqchunkrev)/;
			#}
			#if (defined $bad1) {
			#	my $len = length($bad1);
			#	$seq =~ s/^.{$len}//;
			#	$qua =~ s/^.{$len}//;
			#}
			#if (defined $bad2) {
			#	my $len = length($bad2);
			#	$seq =~ s/.{$len}$//;
			#	$qua =~ s/.{$len}$//;
			#}
			#$bad1 = "" if not defined $bad1;
			#$bad2 = "" if not defined $bad2;
			return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
		}
		else {
			my $bcseq;
			($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) = get_second_barcode($seq, $bc1, $bcseq, $name, $bad, $data, $bc02);
			$barcode3 = $barcode;
			undef $barcode;
		}
	}
	print $outLog "$LRD$name\tNO BARCODE!$N\n" if not defined $barcode;
	$seq_nobc ++ if not defined($barcode);
	return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
}
=cut

__END__


#___
				my $bcseqrev = $bcseq; 
				$bcseqrev =~ tr/ACTG/TGAC/;
				$bcseqrev = rev2($bcseqrev);
				
				print $outLog "\n${YW}1. TEST FORWARD $bcseq$N\n";
				my $muscleinput = ">$bc\n$bcseq\n>$nameshort\n$seq";
				my @res = muscle($muscleinput);
				print "RES:\n$LGN" . join("", @res) . "\n$N\n";
				print "muscleinput:$LPR$muscleinput$N\n";
				($matchperc, $pos) = parse_muscle(\@res, $bcseq, $seq);
				my $rez2 = "";
				for (my $r = 0; $r < @res; $r++) {
					chomp($res[$r]);
					if ($res[$r] =~ /^>/) {
						$res[$r] .= "MUHSPACE";
					}
					$rez2 .= $res[$r];
				}
				$rez2 =~ s/>/\n/g;
				$rez2 =~ s/^\n$//g;
				$rez2 =~ s/MUHSPACE/\n/g;




#			if ($bcseq !~ /,/) {
#				my ($matchperc, $pos, $rez2, $bcseq, $bcseqrev) = testbc($bc, $bcseqorig, $nameshort, $seq);
#			}
#			else {
#				my $seq1 = $seq if length($seq) <= 100;
#				($seq1) = $seq =~ /^(.{100})/ if length($seq) > 100;
#				my ($matchperc, $pos, $rez2, $bcseq, $bcseqrev) = testbc($bc, $bcseqorig, $nameshort, $seq);
#				$seq1 = $seq if length($seq) <= 100;
#				($seq1) = $seq =~ /^(.{100})/ if length($seq) > 100;
#				($matchperc, $pos, $rez2, $bcseq, $bcseqrev) = testbc($bc, $bcseqorig, $nameshort, $seq);
#			}

			#else {
			#	($bctest1, $bctest2) = split(",", $bcseqorig);
			#	$bcseq = $bctest1;
			#}
			#my @res2 = split("\n", $rez2);
			#($pos) = $res2[3] =~ /^(.*)($bcseq|$bcseqrev)/;
			#$pos = defined $pos ? length($pos) : -1;
