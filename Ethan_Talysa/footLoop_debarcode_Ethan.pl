#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use File::Basename qw(dirname); use Cwd qw(abs_path);
use warnings FATAL => 'all';
use vars qw($opt_v $opt_b $opt_i $opt_d $opt_x $opt_2);
getopts("vb:i:d:x2:");

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
my ($barcodeFile2) = $opt_2 if defined $opt_2;
my $bc02;
my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N ${LGN}[-x: Dry run]$N $CY-i <input1>$N $LPR-b <barcode.csv>$N

Optional: -d <name for debug e.g. 203>
";
#-2 <2nd barcode, 2 columns of ID and barcode sequence>

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
open (my $outExtra, ">", "$dbfolder/$inputName\_ExtraBC.tsv") or die "Cannot write to $dbfolder/$inputName\_ExtraBC.tsv: $!\n";
# parse barcode
my ($bc0, %out);
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
	die "\n${LRD}FATAL ERROR$N: Multiple barcode name for the same barcode sequence.\n$LCY" . $bc0->{uc($bcseq)}{bc} . "$N\ncurrent=$LPR$bc$N\nbarcode=$LGN" . uc($bcseq) . "$N\n\n" if defined $bc0->{uc($bcseq)};
	if (@arr == 2) {
		$bc0->{uc($bcseq)}{plasmid} = $plasmid;
		$bc0->{uc($bcseq)}{target} = $target;
		$bc0->{uc($bcseq)}{total} = 0;
		$bc0->{uc($bcseq)}{line} = $linecount;
		$bc0->{uc($bcseq)}{bc} = $bc;
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc,bcseq=$bcseq\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	else {
		$bc0->{uc($bcseq)}{plasmid} = $plasmid;
		$bc0->{uc($bcseq)}{target} = $target;
		$bc0->{uc($bcseq)}{bc} = $bc;
		$bc0->{uc($bcseq)}{total} = 0;
		$bc0->{uc($bcseq)}{line} = $linecount;
		LOG($outLog, "\t$linecount\t$dbfolder/$inputName\_$bc.fastq\tbc=$bcbefore,plasmid=$plasmid,desc=$desc,bcseq=$bcseq\n","NA");
		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	#open ($out{$bc}, ">", "$dbfolder/$inputName\_$bc.fastq") or die "Cannot write to $dbfolder/$inputName\_$bc.fastq: $!\n";
}
print $LGN . "\t" . (1+scalar(@line)) . "\t$dbfolder/UNKNOWN.fastq$N (no barcode)\n\n";

if (defined $opt_2) {
	my @line2 = `cat $barcodeFile2`; $linecount = -1;
	foreach my $line2 (@line2[0..@line2-1]) {
		chomp($line2);
		$linecount ++;
		next if $line2 =~ /^\n$/;
		next if $line2 !~ /[a-zA-Z0-9]+/;
		my ($target2, $bcseq2, $plasmid2) = split("\t", $line2);
		$plasmid2 = uc($plasmid2);
		print $outLog "seq=$bcseq2,plasmid=$plasmid2,desc2=$target2\n";
		$bc02->{uc($bcseq2)}{plasmid} = uc($plasmid2);
		$bc02->{uc($bcseq2)}{target} = $target2;
		$bc02->{uc($bcseq2)}{bc} = uc("$target2");
		$bc02->{uc($bcseq2)}{total} = 0;
		$bc02->{uc($bcseq2)}{line} = $linecount;
		if ($plasmid2 eq "TAC") {
			my $bcseq2revcomp = revcomp($bcseq2);
			print $outLog "bc=$target2 $LGN$bcseq2 $bcseq2revcomp$N\n";
		}
		foreach my $bcseq (sort keys %{$bc0}) {
			my $plasmid = $bc0->{$bcseq}{plasmid};
			#print $outLog "  plasmid=$LGN$plasmid$N\n";
#			if ($plasmid eq "TAC") {print $outLog "bc=$bc0->{$bcseq}{bc} seq=$LGN$bcseq$N\n";}
			next if $plasmid ne $plasmid2;
			if (defined $bc0->{$bcseq}{bc2}) {
				if (defined $bc0->{$bcseq}{bc2}{uc($bcseq2)}) {
					#if ($plasmid eq "TAC") {print $outLog "bc=$bc0->{$bcseq}{bc}\n";}
					die "
Already defined bcseq2=$LCY$bcseq2$N for 
bc=$bc0->{$bcseq}{bc} (bcseq1=$bcseq)
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



if (defined $opt_x) {
	print "\nDry run ${LGN}SUCCESS$N\n\n";
	exit;
}

open (my $outnobc, ">", "$dbfolder/$inputName\_UNKNOWN.fastq") or die "Cannot write to $dbfolder/$inputName\_UNKNOWN.fastq: $!\n";
my ($seq_nobc, $seq_dblbc, $seq_totbc, $readCount) = (0,0,0,0);
my ($data, $bad, @len, $in1);
my ($folder1, $fileName1) = getFilename($input1, "folder");
print "${LGN}1. -b $LCY$barcodeFile$N was parsed successfully!\n$N\n";
# parse fastq
my $wasoptd = 0;
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
#	print "readName=$readName\nreadWant=$opt_d\n";
	if (defined $opt_d) {
		if ($readName eq $opt_d and $wasoptd eq 0) {
			$wasoptd = 1;
		}
		elsif ($readName ne $opt_d and $wasoptd eq 0) {
			next;
		}
		elsif ($readName ne $opt_d and $wasoptd eq 1) {
			last;
		}
	}
	my ($bc, $bc2, $seq2, $qua2, $bc3);
#	return ($barcode, $barcode2, $bc0, $seq, $qua, $data, $barcode3);

	#return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
	print $outLog "\n$YW INFERRING BARCODE$N\n$readCount: $LGN$readName$N\n";
	($bc, $bc2, $bc0, $seq2, $qua2, $data, $bc3, $bc02) = infer_barcode($bc0, $seq, $qua, $readName, $data, $bc02);
	$bc3 = "UNKNOWN" if not defined $bc3;
	$bc = "UNKNOWN" if not defined $bc;
	$bc2 = "UNKNOWN" if not defined $bc2;
	#$bc = "UNKNOWN" if defined $bc and $bc eq "UNDEF";
	print $outLog "$readCount: $LGN$readName$N, bc=$bc, bc3=$bc3\n";
	#print "$readCount: $LGN$readName$N, bc=$bc, bc3=$bc3\n";
	print $outnobc "$readName\n$seq\n\+\n$qua\n" if not defined $seq2;
	print $outLog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
	print $outExtra "$readName\t$bc\t$bc2\t$bc3\n";
	if (defined $seq2) {
		if (not defined $out{$bc}) {
			print "Printed bc=$bc, read=$readName\n";
			open ($out{$bc}, ">", "$dbfolder/$inputName\_$bc3.fastq") or die "Cannot write to $dbfolder/$inputName\_$bc3.fastq: $!\n";
		}
		print {$out{$bc}} "$readName\n$seq2\n\+\n$qua2\n";
	}

	print date() . "\t$YW$input1$N: Done $LGN$readCount$N\n" if $readCount % 500 == 0;
	last if $readCount > 500;
}
close $in1;
close $outExtra;
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
die;
foreach my $bcseq (sort {$bc0->{$a}{line} <=> $bc0->{$b}{line}} keys %{$bc0}) {
	print $outLog "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
	print STDERR "\t$bc0->{$bcseq}{bc}\t$bcseq\t$bc0->{$bcseq}{total}\n";
}

sub get_second_barcode {
	my ($seq, $bc1, $bcseq, $name, $bad, $data, $bc02) = @_;
	print $outLog "\nChecking 2nd barcode\n";
	my $bc2;
	if (not defined $bcseq) {
		$bc2 = $bc02;
	}
	else {
		$bc2 = $bc1->{$bcseq}{bc2};
	}
	my $barcode; my $barcode2; my $barcode3; my $type;
	my $seq1; my $seq2;
	my $seq1chunk; my $seq2chunk;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	for (my $l = 0; $l < 2; $l++) {
		foreach my $bcseq2 (sort keys %{$bc2}) {
			my $bcseq2chunk = substr($bcseq2, $l, length($bcseq2)-$l+1);
			my $length = length($bcseq2chunk);
			my $bcseq2chunkrev = revcomp($bcseq2chunk);
			if ($seq =~ /($bcseq2chunk|$bcseq2chunkrev)/) {
				my $target2 = $bc2->{$bcseq2}{target};
				my $pos = 0 if $seq =~ /^($bcseq2chunk|$bcseq2chunkrev)/;
				   $pos = $seq_length if $seq =~ /($bcseq2chunk|$bcseq2chunkrev)$/;
				my ($pos1, $pos2) = $seq =~ /^(.+)($bcseq2chunk|$bcseq2chunkrev)(.+)$/ if not defined $pos and $seq =~ /($bcseq2chunk|$bcseq2chunkrev)/;
					$pos = defined $pos ? $pos : length($pos1) <= 10 ? length($pos1) : length($pos2);
				my $type = $pos < 0.5 * $seq_length ? 1 : 2;
				my $lengthpos1 = defined $pos1 ? length($pos1) : 0;
				my $lengthpos2 = defined $pos2 ? length($pos2) : 0;
				$pos1 = "POS1_UNDEF" if not defined $pos1; 
				$pos2 = "POS2_UNDEF" if not defined $pos2;
				die "$name: barcode=$barcode, bcseq=$bcseq2 ($bcseq2chunk), pos1=$pos1 pos2=$pos2\n" if not defined $pos;
				print $outLog "get_2nd_barcode: pos = $pos, pos1=" . $lengthpos1 . ", pos2=" . $lengthpos2 . "\nseq: $seq\nbcseq2=$bcseq2\nbcseq2chunk=$bcseq2chunk\n";
				$barcode2 = $seq =~ /$bcseq2chunk/ ? "$type,$bcseq2chunk" : "$type,$bcseq2chunkrev";
				if (defined $barcode) {

					$bad->{$name} = 1;
					$seq2 = $bcseq2;
					$seq2chunk = $bcseq2chunk;
					print $outLog "\tB,l=$l,$LCY$name\tDBL\t$pos/$seq_length\tprev=$barcode\tcurr=$bc2->{$bcseq2}{target}\tprev=$seq1 (chunk=$seq1chunk)\tcurr=$seq2 (chunk=$seq2chunk)$N\n$seq\n";
					#print "\n\n2nd: l=$l, $LCY$name\tDBL\t$pos/$seq_length\tprev=$barcode\tcurr=$bc2->{$bcseq2}{target}\tprev=$seq1 (chunk=$seq1chunk)\tcurr=$seq2 (chunk=$seq2chunk)$N\n$seq\n\n";
				}
				else {
					if (not defined $bcseq) {
						$barcode = $bc2->{$bcseq2}{target};
						$barcode3 = "UNKNOWN";# $barcode;
					}
					else {
						$barcode = $bc1->{$bcseq}{bc};
						$barcode3 = $barcode . "_" . $bc2->{$bcseq2}{target};
						$bc1->{$bcseq}{total} ++;
						$data->{$name}{bc} = $barcode;
						$data->{$name}{pos} = $pos;
					}
#					$barcode3 = $bc2->{$bcseq2}{target};
#					$barcode3 = $barcode . "_" . $bc2->{$bcseq2}{target};
					$seq1 = $bcseq2;
					$seq1chunk = $bcseq2chunk;
					print $outLog "\tB,l=$l,$LGN$name\tFIRST\t$pos/$seq_length\t$barcode$N\n";
					#print "\n\n2nd: l=$l, $LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n\n";
				}
			}
		}
		return ($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) if defined $barcode;
	}
	my $barcodetemp;

	#if (not defined $barcode) {
	($barcodetemp) = $seq =~ /TAGCGG(\w{1,10})AGCACTAA/;
	$barcodetemp .= "_FW" if (defined $barcodetemp);
	if (not defined $barcodetemp) {
		($barcodetemp) = $seq =~ /TTAGTGCT(\w{1,10})CCGCTA/;
		$barcodetemp .= "_RV" if (defined $barcodetemp);
	}

	if (defined $barcodetemp) {
		if (defined $bcseq) {
			$barcode = $bc1->{$bcseq}{bc};
			$barcode2 = $barcodetemp;
			$barcode3 = $barcode . "_" . "UNKNOWN";
		}
		else {
			#$barcode = "UNKNOWN";
			$barcode2 = $barcodetemp;
			$barcode3 = "UNKNOWN"
		}
	}

	return ($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02);
}
sub infer_barcode {
	my ($bc1, $seq, $qua, $name, $data, $bc02) = @_;
	my $barcode; my $barcode2; my $barcode3; my $type;
	my $seq1; my $seq2;
	my $seq1chunk; my $seq2chunk;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	for (my $l = 0; $l < 4; $l++) {
		foreach my $bcseq (sort {$bc1->{$a}{bc} cmp $bc1->{$b}{bc}} keys %{$bc1}) {
			my $bcseqchunk = substr($bcseq, $l, length($bcseq)-$l+1);
			my $length = length($bcseqchunk);
			my $bcseqchunkrev = $bcseqchunk;
			$bcseqchunkrev =~ tr/ACGT/TGCA/;
			$bcseqchunkrev = reverse($bcseqchunkrev);
			#if (($seq =~ /^(.){0,10}($bcseqchunk|$bcseqchunkrev)/ or $seq =~ /($bcseqchunk|$bcseqchunkrev).{0,10}$/)) {
			if ($seq =~ /($bcseqchunk|$bcseqchunkrev)/) {# or $seq =~ /($bcseqchunk|$bcseqchunkrev)/)) {
				if (defined $bc1->{$bcseq}{bc2}) {
					my $typez = $seq =~ /$bcseqchunk/ ? "FW" : "RV";
					print $outLog "\nbcseq=$LGN$bcseq$N, bc=$bc1->{$bcseq}{bc}$N, type=$LPR$typez$N, checking 2nd barcode now\n";
					($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) = get_second_barcode($seq, $bc1, $bcseq, $name, $bad, $data, $bc02);
				}
				if (not defined $barcode) {
					print $outLog "bcseq=$LGN$bcseq has no 2nd barcode!\n";
					my $pos = 0 if $seq =~ /^($bcseqchunk|$bcseqchunkrev)/;
					   $pos = $seq_length if $seq =~ /($bcseqchunk|$bcseqchunkrev)$/;
					my ($pos1, $pos2) = $seq =~ /^(.+)($bcseqchunk|$bcseqchunkrev)(.+)$/ if not defined $pos and $seq =~ /($bcseqchunk|$bcseqchunkrev)/;
						$pos = defined $pos ? $pos : length($pos1) <= 10 ? length($pos1) : length($pos2);
					my $type = $pos < 0.5 * $seq_length ? 1 : 2;
					die "$name: barcode=$barcode, bcseq=$bcseq ($bcseqchunk), pos1=$pos1 pos2=$pos2\n" if not defined $pos;
					#$barcode2 = $seq =~ /$bcseqchunk/ ? "$type,$bcseqchunk" : "$type,$bcseqchunkrev";
					my $lengthpos1 = defined $pos1 ? length($pos1) : 0;
					my $lengthpos2 = defined $pos2 ? length($pos2) : 0;
					$pos1 = "POS1_UNDEF" if not defined $pos1; 
					$pos2 = "POS2_UNDEF" if not defined $pos2;
					print $outLog "infer_barcode: pos = $pos, pos1=" . $lengthpos1 . ", pos2=" . $lengthpos2 . "\nseq: $seq\nbcseq=$bcseq\nbcseqchunk=$bcseqchunk\n";
					if (defined $barcode) {
						$barcode2 = "UNKNOWN" if not defined $barcode2;
						$barcode3 = $barcode . "_" . $barcode2;
						$bad->{$name} = 1;
						$seq2 = $bcseq;
						$seq2chunk = $bcseqchunk;
						print $outLog "\tA,l=$l,$LCY$name\tDBL\t$pos/$seq_length\tprev=$barcode\tcurr=$bc1->{$bcseq}{bc}\tprev=$seq1\tcurr=$seq2$N\n";
						#print "\n\n$LCY$name\tDBL\t$pos/$seq_length\tprev=$barcode\tcurr=$bc1->{$bcseq}{bc}\tprev=$seq1\tcurr=$seq2$N\n\n\n";
					}
					else {
						$barcode = $bc1->{$bcseq}{bc};
						$barcode2 = "UNKNOWN" if not defined $barcode2;
						$barcode3 = $barcode . "_" . $barcode2;
						$seq1 = $bcseq;
						$seq1chunk = $bcseqchunk;
						$data->{$name}{bc} = $barcode3;
						$data->{$name}{pos} = $pos;
						$bc1->{$bcseq}{total} ++;
						print $outLog "\tA,l=$l,$LGN$name\tFIRST\t$pos/$seq_length\t$barcode\t$barcode2\t$barcode3$N\n";
						#print "\n\n$LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n\n\n";
					}
				}
			}
			die if $l > 3;
		}
		$seq_totbc ++ if $l == 0;
		if (defined $barcode and not defined $bad->{$name}) {
			$barcode2 = "UNKNOWN" if not defined $barcode2;
			$barcode3 = $barcode . "_" . $barcode2 if not defined $barcode3;
			$data->{$name}{pos} = "NA" if not defined $data->{$name}{pos};
			#print "$data->{$name}{pos}\n";
			print $outLog "$LGN$name\tGOOD\t$data->{$name}{pos}/$seq_length\t$barcode\t$barcode2\t$barcode3$N\n\n";
			$barcode3 = $barcode if $barcode =~ /pFC9/i;
			$barcode2 = "pFC9" if $barcode =~ /pFC9/i;
			return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
#			my $bcseq = $seq1;
#			my $bcseqchunk    = defined $seq1chunk ? $seq1chu
#			my $bcseqchunkrev = revcomp($seq1chunk);
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
#			return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
		}
#		else {
#			my $bcseq;
#			($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) = get_second_barcode($seq, $bc1, $bcseq, $name, $bad, $data, $bc02);
#			my $barcodeprint = $barcode; $barcodeprint = "UNKNOWN" if not defined $barcode;
#			$data->{$name}{pos} = "NA" if not defined $data->{$name}{pos};
#			$barcode2 = "UNKNOWN" if not defined $barcode2;
#			$barcode3 = $barcodeprint . "_" . $barcode2 if not defined $barcode3;
#			#$barcode3 = $barcode if $barcode =~ /pFC9/i;
#			#$barcode2 = "pFC9" if $barcode =~ /pFC9/i;
#			print $outLog "$LGN$name\tNO BARCODE1\t$data->{$name}{pos}/$seq_length\t$barcodeprint\t$barcode2\t$barcode3$N\n";
#		}
	}
	my $bcseq;
	($barcode, $barcode2, $barcode3, $bad, $bc1, $data, $seq1, $seq1chunk, $bc02) = get_second_barcode($seq, $bc1, $bcseq, $name, $bad, $data, $bc02);
	my $barcodeprint = $barcode; $barcodeprint = "UNKNOWN" if not defined $barcode;
	$data->{$name}{pos} = "NA" if not defined $data->{$name}{pos};
	$barcode2 = "UNKNOWN" if not defined $barcode2;
	$barcode3 = $barcodeprint . "_" . $barcode2 if not defined $barcode3;
	#$barcode3 = $barcode if $barcode =~ /pFC9/i;
	#$barcode2 = "pFC9" if $barcode =~ /pFC9/i;
	print $outLog "$LGN$name\tNO BARCODE1\t$data->{$name}{pos}/$seq_length\t$barcodeprint\t$barcode2\t$barcode3$N\n";
	#print $outLog "$LRD$name\tNO BARCODE!$N\n" if not defined $barcode;
	$seq_nobc ++ if not defined($barcode);
	return ($barcode, $barcode2, $bc1, $seq, $qua, $data, $barcode3, $bc02);
}

__END__
