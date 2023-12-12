#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_i $opt_b $opt_o);
getopts("vi:b:o:");

my ($fastaFile, $barcodeFile, $outputFile) = ($opt_i, $opt_b, $opt_o);
die "\nUsage: $YW$0$N -i $LCY<fastaFile.fa>$N -b $LGN<3 col tsv barcode file>$N [optional: -o $LGN<outputFile>$N]\n\n" unless defined $fastaFile and defined $barcodeFile and -e $fastaFile and -e $barcodeFile;

my ($pos1, $pos2, $pos3) = (20, 8, 23); #primer, .{20}.{8}.{23}, pos2=barcode
my ($neg1, $neg2, $neg3) = (24, 8, 20); #primer, .{20}.{8}.{23}, pos2=barcode
my %data;
my ($fastafolder  , $fastaFilename)   = mitochy::getFilename($fastaFile, "folderfull");
my ($barcodefolder, $barcodeFilename) = mitochy::getFilename($barcodeFile, "folderfull");
$outputFile = "$fastaFilename\_$barcodeFilename\_bigfasta.fa" if not defined $outputFile;

my $linecount = 0;
my $inputFA = $fastaFile;
open (my $inFASTA, "<", $inputFA) or die "Cannot read from $inputFA: $!\n";
my $fasta = new FAlite($inFASTA);
while (my $entry = $fasta->nextEntry()) {
	$linecount ++;
	my $seq = $entry->seq;
	my $def = $entry->def;
	$def =~ s/^>//;
	$def = uc($def);
	$data{$def}{seq} = $seq;
}
close $inFASTA;
my %def;
open (my $out1, ">", $outputFile) or die "Cannot write to $outputFile: $!\n";
open (my $out2, ">", "$outputFile.tsv") or die "Cannot write to $outputFile.tsv: $!\n";
open (my $in1, "<", $barcodeFile) or die "Failed to read from $barcodeFile: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if $line =~ /^#/;
	my ($id, $plasmid, $desc, $posseq, $negseq) = @arr;
	$plasmid = uc($plasmid);
	my ($pos1seq, $pos2seq, $pos3seq) = $posseq =~ /^(.{$pos1})(.{$pos2})(.{$pos3})$/;
	my ($neg1seq, $neg2seq, $neg3seq) = $negseq =~ /^(.{$neg1})(.{$neg2})(.{$neg3})$/;
	my $count = 0;
	if ($plasmid =~ /mix/i) {
		my ($plasmid2) = $plasmid =~ /^(.+)mix$/i;
		foreach my $def (sort keys %data) {
			next if $def =~ /ApaLI/i;
			next if $def !~ /^$plasmid2/;
			$count ++;
			my $seq = $data{$def}{seq};
			die "pos: Seq does not contain expected primers!\n$seq\n${pos1seq}NNNNNNNN${pos3seq}\n$pos1seq$pos2seq$pos3seq\n\n" if $seq !~ /${pos1seq}NNNNNNNN${pos3seq}/;
			die "neg: Seq does not contain expected primers!\n$seq\n${neg1seq}NNNNNNNN${neg3seq}\n$neg1seq$neg2seq$pos3seq\n\n" if $seq !~ /${neg1seq}NNNNNNNN${neg3seq}/;
			$seq =~ s/${pos1seq}NNNNNNNN${pos3seq}/${pos1seq}${pos2seq}${pos3seq}/;
			$seq =~ s/${neg1seq}NNNNNNNN${neg3seq}/${neg1seq}${neg2seq}${neg3seq}/;
			#print ">id$id\_def$def\_plasmid$plasmid\_fwprimer$posseq\_rvprimer$negseq\_desc$desc\_fwbc$pos2seq\_rvbc$neg2seq\n$seq\n";
			my $newdef = "$def\_desc$desc";
			my $newdefcount = 1;
			while (1) {
				last if not defined $def{$newdef};
				$newdefcount ++;
				$newdef = "$def\_desc$desc\_rep$newdefcount";
			}
			$def{$newdef} = 1;
			print $out2 "$id\t$def\t$plasmid\t$posseq\t$negseq\t$desc\t$pos2seq\t$neg2seq\t$newdef\t$seq\n";
			print $out1 ">$newdef\n$seq\n";
#			print "\t$def\n";
		}
	}
	else {
		$count ++ if defined $data{$plasmid};
		my $seq = $data{$plasmid}{seq};
		die "pos: Seq does not contain expected primers!\n$seq\n${pos1seq}NNNNNNNN${pos3seq}\n$pos1seq$pos2seq$pos3seq\n\n" if $seq !~ /${pos1seq}NNNNNNNN${pos3seq}/;
		die "neg: Seq does not contain expected primers!\n$seq\n${neg1seq}NNNNNNNN${neg3seq}\n$neg1seq$neg2seq$pos3seq\n\n" if $seq !~ /${neg1seq}NNNNNNNN${neg3seq}/;
		$seq =~ s/${pos1seq}NNNNNNNN${pos3seq}/${pos1seq}${pos2seq}${pos3seq}/;
		$seq =~ s/${neg1seq}NNNNNNNN${neg3seq}/${neg1seq}${neg2seq}${neg3seq}/;

		my $newdef = "$plasmid\_desc$desc";
		my $newdefcount = 1;
		while (1) {
			last if not defined $def{$newdef};
			$newdefcount ++;
			$newdef = "$plasmid\_desc$desc\_rep$newdefcount";
		}
		$def{$newdef} = 1;

		print $out2 "$id\t$plasmid\t$plasmid\t$posseq\t$negseq\t$desc\t$pos2seq\t$neg2seq\t$newdef\t$seq\n";
		print $out1 ">$newdef\n$seq\n";
#			print "\t$def\n";
#		print "$id\t$plasmid\t$plasmid\t$posseq\t$negseq\t$desc\t$pos2seq\t$neg2seq\n";
	}
	die "Undefined $plasmid\n" if $count eq 0;
}
close $out1;
close $out2;
