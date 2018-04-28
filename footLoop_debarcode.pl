#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_b $opt_i $opt_d);
getopts("vb:i:d:");
my $usage = "\nUsage: $YW$0$N $CY-i <input1>$N $LPR-b <barcode.csv>$N\nOptional: -d <name for debug e.g. 203>\n";
my $debugname = $opt_d;
$debugname = "NA" if not defined $opt_d;
my ($barcodeFile, $input1) = ($opt_b, $opt_i);
print $usage if not defined $input1 or not defined $barcodeFile;
print "$LRD ERROR$N -i not defined\n" if not defined $input1;
print "$LRD ERROR$N -b not defined\n" if not defined $barcodeFile;
die "$usage\n\nDoes not exist -b ${LPR}barcodeFile$N!\n" if defined $barcodeFile and not -e $barcodeFile;
print "-b: $LPR$barcodeFile$N\n" if defined $barcodeFile;
die "$usage\n\nDoes not exist -i $LCY$input1$N!\n" if defined $input1 and not -e $input1;
print "-i: $LCY$input1$N\n" if defined $input1;
die if not defined $input1 or not defined $barcodeFile or (defined $input1 and not -e $input1) or (defined $barcodeFile and not -e $barcodeFile);

my ($inputName) = getFilename($input1);
my $dbfolder = "debarcode_result";
mkdir $dbfolder if not -d $dbfolder;

open (my $outlog, ">", "$dbfolder/$inputName\_logfile.txt") or die "Cannot write to $dbfolder/logfile.txt: $!\n";
# parse barcode
my (%bc, %out);
my @line = `cat $barcodeFile`; my $linecount = 0;
print "\nOutput:\n";
foreach my $line (@line[0..@line-1]) {
	chomp($line);
	$linecount ++;
	next if $line =~ /^\n$/;
	next if $line !~ /[a-zA-Z0-9]+/;
	next if $line =~ /^No/ or $line =~ /^target/i;
	my $delim = $line =~ /\t/ ? "\\t" : $line =~ /,/ ? "," : $line =~ /[ ]/ ? " " : die "Cannot determine delimiter at line=$LCY\n$line\n$N\n";
	my @arr = split("$delim", $line);
	my ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;
	if (@arr == 2) {
		($bc, $bcseq) = @arr;
		
	}
	else {
		($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#split(",", $line)  if $line !~ /\t/;
	   ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;#osplit("\t", $line) if $line =~ /\t/;
		$desc =~ s/[ ]+//g;
		$plasmid =~ s/[ ]+//g;
	}
	$bcseq =~ s/^[ \s]+//g;
	$bcseq =~ s/[ \s]+$//g;
	$bc =~ s/^[ \s]+//g;
	$bc =~ s/[ \s]+$//g;
	$bc = "$bc\_$plasmid\_$desc" if defined $plasmid and defined $desc;
	$plasmid = "" if not defined $plasmid;
	$desc = "" if not defined $desc;
	$bc =~ s/_$//g;
	$bc =~ s/^_//g;
	print "$bc\t$bcseq\n";
	$target = $bc if not defined $target;
	die "Multiple barcode name for the same barcode sequence (" . uc($bcseq) . ")\n" if defined $bc{uc($bcseq)};
	if (@arr == 2) {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		$bc{uc($bcseq)}{bc} = $bc;
	}
	else {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{bc} = $bc;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		print $outlog "$linecount\t$dbfolder/$inputName\_$bc.fastq\n";
		print STDERR "$LGN$linecount\t$dbfolder/$inputName\_$bc.fastq$N\n";
	}
	open ($out{$bc}, ">", "$dbfolder/$inputName\_$bc.fastq") or die "Cannot write to $dbfolder/$inputName\_$bc.fastq: $!\n";
}
print $LGN . scalar(@line) . "\t$dbfolder/UNKNOWN.fastq$N (no barcode)\n\n";
open (my $outnobc, ">", "$dbfolder/$inputName\_UNKNOWN.fastq") or die "Cannot write to $dbfolder/$inputName\_UNKNOWN.fastq: $!\n";
my ($seq_nobc, $seq_dblbc, $seq_totbc, $readCount) = (0,0,0,0);
my (%data, %bad, @len, $in1);
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
print "${LGN}1. -b $LCY$barcodeFile$N was parsed successfully!\n$N\n";
# parse fastq
my $isDir = -d $input1 ? "(is a directory!)" : -e $input1 ? "" : "($input1 does not exist!)";
die "$LRD FATAL ERROR$N: -i $LCY$input1$N is not a fastq file! $LRD$isDir$N\n" if -d $input1 or not -e $input1;
open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.gz$/i;
open ($in1, "zcat $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.gz$/i;
while (my $line = <$in1>) { chomp($line);
	# 1. name: m171023_184836_42145_c101402632550000001823289101251862_s1_p0/0/ccs
	$readCount ++;
	my ($readName) = $line;
#	my ($readName) = $line =~ /\/(\w+)\/ccs/;
#       $readName  = $line if not defined $readName;
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
	print $outlog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
	print {$out{$bc}} "$readName\n$seq2\n\+\n$qua2\n" if defined $seq2;

	print date() . "\t$YW$input1$N: Done $LGN$readCount$N\n" if $readCount % 500 == 0;
#DEBUG	last if $readCount % 1000 == 0;
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
print $outlog "

Total seq : $seq_totbc
Has bc    : $seq_goodbc ($seq_goodbcPerc \%)
Double bc : $double_bc
No bc     : $seq_nobc

";
foreach my $bcseq (sort {$bc{$a}{line} <=> $bc{$b}{line}} keys %bc) {
	print $outlog "\t$bc{$bcseq}{bc}\t$bc{$bcseq}{total}\n";
	print STDERR "\t$bc{$bcseq}{bc}\t$bc{$bcseq}{total}\n";
}
#print $outlog "\n";\n###print STDERR "\n";

sub infer_barcode {
	my ($seq, $qua, $name) = @_;
	my $barcode; my $barcode2; my $type;
	my $seq1; my $seq2;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	for (my $l = 0; $l < 4; $l++) {
#		last if $l != 0;
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
					#print "$LRD DBL$N  $name\t$l\t$bcseq2: $bcseq\n" if $l != 0;
					$bad{$name} = 1;
					$seq2 = $bcseq;
					print $outlog "$LCY$name\tDBL\t$pos/$seq_length\t$barcode\t$bc{$bcseq2}{bc}\t$seq1\t$seq2$N\n";
				}
				else {
					#print "$LGN GOOD$N $name\t$l\t$bcseq2: $bcseq\n" if $l != 0;
					$barcode = $bc{$bcseq2}{bc};
					$seq1 = $bcseq;
					$data{$name} = $barcode;
					$bc{$bcseq2}{total} ++;
					print $outlog "$LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n";
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
#			print "$name, barcode = $LGN$bcseq$N or $LCY$bcseqrev$N\n$seq\n$qua\n";
			if (defined $bad1) {
				my $len = length($bad1);
#				print "-> beginning has bad of length $LGN$len$N\n";
				$seq =~ s/^.{$len}//;
				$qua =~ s/^.{$len}//;
			}
			if (defined $bad2) {
				my $len = length($bad2);
#				print "-> end has bad of length $LCY$len$N\n";
				$seq =~ s/.{$len}$//;
				$qua =~ s/.{$len}$//;
			}
			$bad1 = "" if not defined $bad1;
			$bad2 = "" if not defined $bad2;
#			print "\n$LGN$bad1$N$seq$LCY$bad2$N\n$LGN$bad1$N$qua$LCY$bad2$N\n";
			return ($barcode, $barcode2, $seq, $qua);
		}
	}
	print $outlog "$LRD$name\tNO BARCODE!$N\n" if not defined $barcode;
	$seq_nobc ++ if not defined($barcode);
	return;
}

__END__
==line 198
	for (my $l = 1; $l <= 3; $l++) {
	foreach my $bcseq2 (sort keys %bc) {
		my $bcseq = substr($bcseq2, $l, length($bcseq2)-$l+1);
		my $length = length($bcseq);
		my $bcseqrev = $bcseq;
		$bcseqrev =~ tr/ACGT/TGCA/;
		$bcseqrev = reverse($bcseqrev);
		if (not defined $barcode and ($seq =~ /$bcseq/ or $seq =~ /$bcseqrev/)) {
			my $pos = 0 if $seq =~ /^$bcseq/ or $seq =~ /^$bcseqrev/;
			   $pos = $seq_length if $seq =~ /$bcseq$/ or $seq =~ /$bcseqrev$/;
			my ($pos1, $pos2) = $seq =~ /^(.+)$bcseq(.+)$/ if not defined $pos and $seq =~ /$bcseq/;
			   ($pos1, $pos2) = $seq =~ /^(.+)$bcseqrev(.+)$/ if not defined $pos and $seq =~ /$bcseqrev/;
				$pos = defined $pos ? $pos : $seq =~ /$bcseq/ ? length($pos1) : $seq =~ /$bcseqrev/ ? length($pos2) : "-1";
			my $type = $pos < 0.5 * $seq_length ? 1 : 2;
			$barcode2 = $seq =~ /$bcseq/ ? "$type,$bcseq" : "$type,$bcseqrev";
			die "$name: barcode=$barcode, bcseq=$bcseq, pos1=$pos1 pos2=$pos2\n" if not defined $pos;
			$barcode = $bc{$bcseq2}{bc};
			$seq1 = $bcseq;
			$data{$name} = $barcode;
			$bc{$bcseq2}{total} ++;
			print $outlog "$LCY$name\tOffseq=$l; GOOD\t$pos/$seq_length\t$barcode$N\n";
		}
		elsif (defined $barcode and ($seq =~ /$bcseq/ or $seq =~ /$bcseqrev/) and $barcode ne $bc{$bcseq2}{bc}) {
			my $pos = 0 if $seq =~ /^$bcseq/ or $seq =~ /^$bcseqrev/;
			   $pos = $seq_length if $seq =~ /$bcseq$/ or $seq =~ /$bcseqrev$/;
			my ($pos1, $pos2) = $seq =~ /^(.+)$bcseq(.+)$/ if not defined $pos and $seq =~ /$bcseq/;
			   ($pos1, $pos2) = $seq =~ /^(.+)$bcseqrev(.+)$/ if not defined $pos and $seq =~ /$bcseqrev/;
				$pos = defined $pos ? $pos : $seq =~ /$bcseq/ ? length($pos1) : $seq =~ /$bcseqrev/ ? length($pos2) : "-1";
			my $type = $pos < 0.5 * $seq_length ? 1 : 2;
			$barcode2 .= $seq =~ /$bcseq/ ? ";$type,$bcseq" : ";$type,$bcseqrev";
			$bad{$name} = 1;
			$seq2 = $bcseq if $seq =~ /$bcseq/;
			$seq2 = $bcseqrev if $seq =~ /$bcseqrev/;
			print $outlog "$LPR$name\tOffset=$l; DBL\t$pos/$seq_length\t$barcode\t$bc{$bcseq2}{bc}\t$seq1\t$seq2$N\n";
		}
	}
	last if defined $barcode;
	}
