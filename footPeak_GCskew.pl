#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_i);
getopts("vi:");

die "\nusage: $YW$0$N -i $CY<input1>$N\n\n" unless defined $opt_i and -d $opt_i;
my ($folder) = $opt_i;

my %gene;
my @line = `cat /group/stella/Work/Data/Fastq/110101_pacbio/1_Bed/180407_task1_GCskew/0_all_indexes_Alpha301.bed`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene, $zero, $strand, $feature) = split("\t", $line);
	die "Undefine gene at line=$line\n" if not defined $gene;
	$gene{$gene}{feature} = $feature;
}

my @files = <$folder/*.tsv>;
my %data;
my @header = ("label", "gene", "strand", "window", "threshold", "convtype", "wind2", "sample", "type");
print "\n\nThere i no file with .tsv in $LCY$folder$N!\n" and exit if (@files == 0);
foreach my $input1 (sort @files) {
	my ($WINDOW, $SAMPLE, $TYPE);
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");
	#my ($label, $barcode, $desc, $gene, $strand, $window, $threshold, $convtype, $wind2, $sample, $type) = $fileName1 =~ /^(PCB.+)_(BC\d+)?_?(\w+)?_?gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	my @arr = $fileName1 =~ /^(PCB\w+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CG|CH|GH|GC).PEAK.genome.bed_(\d+)_([A-Z]).temp.fa.(\w+).tsv/;
	$arr[0] =~ s/^(PCB\d+)_.+$/$1/;
	my $outName = join("_", @arr[0..6]) . "_" . $arr[8];
	for (my $i = 0; $i < @arr; $i++) {
		die "Undefined i=$i header=$header[$i] arr[i] undef\n" if not defined $arr[$i];# and $header[$i] !~ /(barcode|desc)/;
		$data{data}{$outName}{$header[$i]} = $arr[$i];
		$WINDOW = $arr[$i] if $header[$i] eq "wind2";
		$SAMPLE = $arr[$i] if $header[$i] eq "sample";
		$TYPE = $arr[$i] if $header[$i] eq "type";
	}
	print "$input1:$LCY$outName$N\n";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		my ($read, $value) = split("\t", $line);
		$data{read}{$outName}{$read}{$SAMPLE} = $value;
		$data{input}{$outName}{$SAMPLE} = join("_", @arr[0..5]);
	}	
	close $in1;
}

open (my $out1, ">", "RESULT.TSV") or die "Cannot write to RESULT.TSV: $!\n";
foreach my $outName (sort keys %{$data{input}}) {
	#open (my $out1, ">", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	print $out1 "file\tread\twindow\ttype\tfeature";
	foreach my $sample (sort keys %{$data{input}{$outName}}) {
		print $out1 "\t$sample";
	}
	print $out1 "\n";
	last;
	#close $out1;
}

foreach my $outName (sort keys %{$data{read}}) {
#open (my $out1, ">>", "$outName.TSV") or die "Cannot write to $outName.TSV: $!\n";
	my $WINDOW = $data{data}{$outName}{wind2};
	my $TYPE = $data{data}{$outName}{type};
	my $SAMPLE = $data{data}{$outName}{sample};
	my $GENE = $data{data}{$outName}{gene};
	my $feature = $gene{$GENE}{feature}; print "Undef gene=$GENE feature\n" and next if not defined $feature;
	foreach my $read (sort keys %{$data{read}{$outName}}) {
		print $out1 "$data{input}{$outName}{$SAMPLE}\t$read\t$WINDOW\t$TYPE\t$feature";
		foreach my $sample (sort keys %{$data{read}{$outName}{$read}}) {
			print $out1 "\t$data{read}{$outName}{$read}{$sample}";
		}
		print $out1 "\n";
	}
	#close $out1;
}
__END__



close $out1;

__END__
PCB1_geneFUS_Pos_20_0.65_CH.PEAK.genome.bed_100_E.temp.fa.dens.tsv
