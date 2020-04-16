#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use File::Basename qw(dirname); use Cwd qw(abs_path);
use vars qw($opt_v $opt_b $opt_i $opt_d $opt_x);
getopts("vb:i:d:x");

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

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N ${LGN}[-x: Dry run]$N $CY-i <input1>$N $LPR-b <barcode.csv>$N

Optional: -d <name for debug e.g. 203>

";

print $usage  unless defined $opt_i and -e $opt_i and defined $opt_b;
print $usage if not defined $opt_i or not defined $opt_b;
print "$LRD ERROR$N -i not defined\n" if not defined $opt_i;
print "$LRD ERROR$N -b not defined\n" if not defined $opt_b;

die "$usage\n\nDoes not exist -b ${LPR}opt_b$N!\n" if not defined $opt_b;
print "-b: $LPR$opt_b$N\n" if defined $opt_b;
die "$usage\n\nDoes not exist -i $LCY$opt_i$N!\n" if defined $opt_i and not -e $opt_i;
print "-i: $LCY$opt_i$N\n" if defined $opt_i;
die "\n\n" if not defined $opt_x and (not defined $opt_i or not defined $opt_b or (defined $opt_i and not -e $opt_i));
die "\n\n" if defined $opt_x and (not defined $opt_b);

my $input1 = $opt_i;
my $debugname = $opt_d;
$debugname = "NA" if not defined $opt_d;
my $barcodeFile = $opt_b;
my @barcodeFiles = split(",", $barcodeFile);

if (defined $opt_x and not defined $opt_i) {
	system("touch dummy.temp") == 0 or die "Failed to create $LGN dummy.temp$N: $!\n";
	$opt_i = "dummy.temp";
	$input1 = $opt_i;
}

my ($inputName) = getFilename($input1);
my $dbfolder = "debarcode_result";
mkdir $dbfolder if not -d $dbfolder;

open (my $outLog, ">", "$dbfolder/$inputName\_logfile.txt") or die "Cannot write to $dbfolder/logfile.txt: $!\n";
# parse barcode
my (%bc, %out);
my @line = `cat @barcodeFiles`; my $linecount = -1;
print "\nOutput:\n";
foreach my $line (@line[0..@line-1]) {
	chomp($line);
	$linecount ++;
	next if $line =~ /^\n$/;
	next if $line !~ /[a-zA-Z0-9]+/;
	next if $line =~ /^(\#|No)/i or $line =~ /\ttarget/i;
	#my $delim = $line =~ /\t/ ? "\\t" : $line =~ /,/ ? "," : $line =~ /[ ]/ ? " " : die "Cannot determine delimiter at line=$LCY\n$line\n$N\n";
	#print "Delimiter = $delim\n";
	my @arr = split("[\t ]+", $line);
	my ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;
	if (@arr == 2) {
		($bc, $bcseq, $desc, $plasmid) = @arr;
	}
	else {
		($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#split(",", $line)  if $line !~ /\t/;
	   ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#osplit("\t", $line) if $line =~ /\t/;
		$desc =~ s/[ ]+//g;
		$plasmid =~ s/[ ]+//g;
	}
	$desc = "NA" if not defined $desc;
	$plasmid = "NA" if not defined $plasmid;
	$target = $bc if not defined $target;
	$desc =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$bcseq =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$bc =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$plasmid =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$target =~ s/(^[ \s]+|[ \s]+$|[\_\-\=\+\/\.\,\<\>\?\'\;\"\:\]\[\}\{\`\~\!\@\#\$\%\^\&\*\(\)\\]+)//g;
	$bc = uc($bc);
	$plasmid = uc($plasmid);
	$desc = uc($desc);
	my $bcbefore = $bc;
	$bc = "bc$bc\_plasmid$plasmid\_desc$desc" if $plasmid ne "NA" or $desc ne "NA";#defined $plasmid and defined $desc;
	die "\n${LRD}FATAL ERROR$N: Multiple barcode name for the same barcode sequence.\n$LCY" . $bc{uc($bcseq)}{bc} . "$N\ncurrent=$LPR$bc$N\nbarcode=$LGN" . uc($bcseq) . "$N)\n\n" if defined $bc{uc($bcseq)};
	if (@arr == 2) {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		$bc{uc($bcseq)}{bc} = $bc;
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	else {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{bc} = $bc;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	open ($out{$bc}, ">", "$dbfolder/$inputName\_$bc.fastq") or die "Cannot write to $dbfolder/$inputName\_$bc.fastq: $!\n";
}
print $LGN . "\t" . (1+scalar(@line)) . "\t$dbfolder/UNKNOWN.fastq$N (no barcode)\n\n";

if (defined $opt_x) {
	print "\nDry run ${LGN}SUCCESS$N\n\n";
	exit;
}

open (my $outnobc, ">", "$dbfolder/$inputName\_UNKNOWN.fastq") or die "Cannot write to $dbfolder/$inputName\_UNKNOWN.fastq: $!\n";
my ($seq_nobc, $seq_dblbc, $seq_totbc, $readCount) = (0,0,0,0);
my (%data, %bad, @len, $in1);
my ($folder1, $fileName1) = getFilename($input1, "folder");
print "${LGN}1. -b $LCY$barcodeFile$N was parsed successfully!\n$N\n";
# parse fastq
my $isDir = -d $input1 ? "(is a directory!)" : -e $input1 ? "" : "($input1 does not exist!)";
die "$LRD FATAL ERROR$N: -i $LCY$input1$N is not a fastq file! $LRD$isDir$N\n" if -d $input1 or not -e $input1;
open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.gz$/i;
open ($in1, "zcat < $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.gz$/i;
while (my $line = <$in1>) { chomp($line);
	# 1. name: m171023_184836_42145_c101402632550000001823289101251862_s1_p0/0/ccs
	$readCount ++;
	my ($readName) = $line;
	$line = <$in1>; chomp($line);
	# 2. seqs: ACTGTCG....
	my $len = int(length($line));
	my ($seq) = $line;
	
	$line = <$in1>; chomp($line);
	# 3. junk (+)

	$line = <$in1>; chomp($line);
	# 4. Quality
	my $qua = $line;

	my ($bc, $bcs, $seq2, $qua2) = infer_barcode($seq, $qua, $readName);
	print $outnobc "$readName\n$seq\n\+\n$qua\n" if not defined $seq2;
	print $outLog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
	print {$out{$bc}} "$readName\n$seq2\n\+\n$qua2\n" if defined $seq2;

	print date() . "\t$YW$input1$N: Done $LGN$readCount$N\n" if $readCount % 500 == 0;
}
close $in1;
print "${LGN}2. -i $LCY$input1$N was parsed successfully!\n$N\n";
my $double_bc = (keys %bad);
my $seq_goodbc = $seq_totbc - $seq_nobc - $double_bc;
my $seq_goodbcPerc = int($seq_goodbc / $seq_totbc * 1000)/10;
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
foreach my $bcseq (sort {$bc{$a}{line} <=> $bc{$b}{line}} keys %bc) {
	print $outLog "\t$bc{$bcseq}{bc}\t$bcseq\t$bc{$bcseq}{total}\n";
	print STDERR "\t$bc{$bcseq}{bc}\t$bcseq\t$bc{$bcseq}{total}\n";
}

sub infer_barcode {
	my ($seq, $qua, $name) = @_;
	my $barcode; my $barcode2; my $type;
	my $seq1; my $seq2;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	for (my $l = 0; $l < 4; $l++) {
		foreach my $bcseq2 (sort {$bc{$a}{bc} cmp $bc{$b}{bc}} keys %bc) {
			my $bcseq = substr($bcseq2, $l, length($bcseq2)-$l+1);
			my $length = length($bcseq);
			my $bcseqrev = $bcseq;
			$bcseqrev =~ tr/ACGT/TGCA/;
			$bcseqrev = reverse($bcseqrev);
			if (($seq =~ /^(.){0,10}($bcseq|$bcseqrev)/ or $seq =~ /($bcseq|$bcseqrev).{0,10}$/)) {
				my $pos = 0 if $seq =~ /^($bcseq|$bcseqrev)/;
				   $pos = $seq_length if $seq =~ /($bcseq|$bcseqrev)$/;
				my ($pos1, $pos2) = $seq =~ /^(.+)($bcseq|$bcseqrev)(.+)$/ if not defined $pos and $seq =~ /($bcseq|$bcseqrev)/;
					$pos = defined $pos ? $pos : length($pos1) <= 10 ? length($pos1) : length($pos2);
				my $type = $pos < 0.5 * $seq_length ? 1 : 2;
				die "$name: barcode=$barcode, bcseq=$bcseq, pos1=$pos1 pos2=$pos2\n" if not defined $pos;
				$barcode2 = $seq =~ /$bcseq/ ? "$type,$bcseq" : "$type,$bcseqrev";
				if (defined $barcode) {
					$bad{$name} = 1;
					$seq2 = $bcseq;
					print $outLog "$LCY$name\tDBL\t$pos/$seq_length\t$barcode\t$bc{$bcseq2}{bc}\t$seq1\t$seq2$N\n";
				}
				else {
					$barcode = $bc{$bcseq2}{bc};
					$seq1 = $bcseq;
					$data{$name} = $barcode;
					$bc{$bcseq2}{total} ++;
					print $outLog "$LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n";
				}
			}
			die if $l > 3;
		}
		$seq_totbc ++ if $l == 0;
		if (defined $barcode and not defined $bad{$name}) {
			my $bcseq    = $seq1;
			my $bcseqrev = $bcseq;
			$bcseqrev =~ tr/ACGT/TGCA/;
			$bcseqrev = reverse($bcseqrev);
			my ($bad1, $good1) = $seq =~ /^(.{0,14}$bcseq|.{0,14}$bcseqrev)(.+)$/ if $seq =~ /^(.{0,14}$bcseq|.{0,14}$bcseqrev)/;
			my ($good2, $bad2) = $seq =~ /^(.+)($bcseq.{0,14}|$bcseqrev.{0,14})$/ if $seq =~ /^(.+)($bcseq.{0,14}|$bcseqrev.{0,14})$/;
			if (defined $bad1) {
				my $len = length($bad1);
				$seq =~ s/^.{$len}//;
				$qua =~ s/^.{$len}//;
			}
			if (defined $bad2) {
				my $len = length($bad2);
				$seq =~ s/.{$len}$//;
				$qua =~ s/.{$len}$//;
			}
			$bad1 = "" if not defined $bad1;
			$bad2 = "" if not defined $bad2;
			return ($barcode, $barcode2, $seq, $qua);
		}
	}
	print $outLog "$LRD$name\tNO BARCODE!$N\n" if not defined $barcode;
	$seq_nobc ++ if not defined($barcode);
	return;
}

__END__
