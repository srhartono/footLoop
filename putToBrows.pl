#!/usr/bin/perl

use warnings;use strict;use Getopt::Std;
use Cwd qw(abs_path);use File::Basename qw(dirname);
use vars qw($opt_i $opt_f $opt_x $opt_y $opt_o $opt_p $opt_t);
getopts("i:f:x:y:o:p:t:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib;use FAlite;

my $usage = "Usage: -i <BED file name> -f <genome.fa e.g. hg19.fa> -x <[0] beg buffer> -y <[0] end buffer> -o <Output file name> -p <Conversion threshold used> -t <track name>\n";
die $usage if not defined $opt_i or not -e $opt_i or not defined $opt_f or not -e $opt_f or not defined $opt_o or not defined $opt_p or not defined $opt_t;

my ($bufferL) = defined $opt_x ? $opt_x : 0;
my ($bufferR) = defined $opt_y ? $opt_y : 0;
die "Left (-x $bufferL) and Right (-y $bufferR) buffer must be integer!\n" if $bufferL !~ /^\-?\d+$/ or $bufferR !~ /^\-?\d+$/;

my $output = $opt_o;$output =~ s/\/+/\//g;
my @folders = split("/", $output);
$output = "";
print "Output0 = $output\n";
for (my $i = 0;$i < @folders-1;$i++) {
	$output .= $folders[$i] . "/";
	$output =~ s/^\/+/\//;
	if (not -d $output) {
		print localtime() . ": Creating Directory: $output\n";
		system("mkdir $output") == 0 or die "Failed to create directory $CY$output$N: $!\n";
	}
}
#die "Output folder $output does not exist somehow (permission issue?)\n" unless -d $output;
$output .= $folders[@folders-1];


my $GTFPos = $output . "_pos.gtf";
open (my $GTFPosfile, ">", $GTFPos) or die "Cannot write to $GTFPos: $!\n";

my $GTFNeg = $output . "_neg.gtf";
open (my $GTFNegfile, ">", $GTFNeg) or die "Cannot write to $GTFNeg: $!\n";

my $bedFile = $opt_i;

my %data;
my @reads;
my @bedLines;
open (my $bedIn, "<", $bedFile) or die "Cannot read from $bedFile: $!\n";
while (my $bedLine = <$bedIn>) {
	chomp($bedLine);
	my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $bedLine);
	if (not defined $strand) {
		$strand = $zero;
	}
	%{$data{$gene}} = (
		'chr' => $chr,
		'beg' => $beg,
		'end' => $end,
		'strand' => $strand,
		'posFile' => "$gene\_Pos$opt_p.txt",
		'negFile' => "$gene\_Neg$opt_p.txt"
	);
	print "Gene = $gene, $chr:$beg-$end ($strand)\n";
}
close $bedIn;

my $header = parse_fasta($output); my %head = %{$header};

foreach my $gene (sort keys %data) {
	my $posFile = $data{$gene}{posFile};
	my $negFile = $data{$gene}{negFile};
	my $pname = $opt_t . "Pos";
	my $nname = $opt_t . "Neg";
	print $GTFPosfile "track name=$pname color=255,0,0\n";
	print $GTFNegfile "track name=$nname color=0,0,255\n";
	parse($gene, $posFile, $data{$gene}, $GTFPosfile) if -e $posFile and -s $posFile != 0;
	print "\tgene $gene\ 's pos file $posFile does not exist or has nothing in it!\n" if not -e $posFile or (-e $posFile and -s $posFile == 0);
	parse($gene, $negFile, $data{$gene}, $GTFNegfile) if -e $negFile and -s $negFile != 0;
	print "\tgene $gene\'s neg file $negFile does not exist or has nothign in it!\n" if not -e $negFile or (-e $negFile and -s $negFile == 0);
}

sub parse_fasta {
	my ($output) = @_;
	# get fasta
	my $fasta = $output . "fa";my $fastaCount = 0;
	while (-e $fasta) {
		$fastaCount ++;
		$fasta = $output . $fastaCount . "fa";
	}
	#not stranded because we agree to use relative to + strand
	my %head;
	my $facmd = "fastaFromBed -fi $opt_f -bed $opt_i -fo $fasta -name";
	system($facmd) == 0 or die "Failed to run $CY$facmd$N: $!\n";
	open (my $faIn, "<", $fasta) or die "Failed to read from $fasta: $!\n";
	$fasta = new FAlite($faIn);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;$def =~ s/^>//;my $gene = $def;
		my @seq = split("", uc($entry->seq));
		my $beg = $data{$gene}{beg};die if not defined $beg;
		my ($beg0, $beg1) = ($beg, $beg + 1);#($beg - 500, $beg - 500+ 1);$beg0 = 1 if $beg0 < 0;$beg1 = 1 if $beg1 < 1;
		my $end = $data{$gene}{end};die if not defined $end;
		my ($end0, $end1) = ($end - 1, $end); #($end - 500, $end - 500+ 1);$end0 = 1 if $end0 < 0;$end1 = 1 if $end1 < 1;
		my $chr = $data{$gene}{chr};die if not defined $chr;
		my $strand = $data{$gene}{strand};die if not defined $strand;
		$head{BLOCK}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.BLOCK\"; transcript_id \"MYNUMBER.$gene.BLOCK\"\n";
		$head{BLOCK}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.BLOCK\"; transcript_id \"MYNUMBER.$gene.BLOCK\"\n";
		$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$beg0\t$beg1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		for (my $i = 0;$i < @seq;$i++) {
			my $begPeak = $i + $beg;# - $bufferL;# 0 based
			my $endPeak = $begPeak + 1;# 1 based
			if ($seq[$i] eq "C" and defined $seq[$i+1] and $seq[$i+1] eq "G") {
				$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t" . ($begPeak+1) . "\t" . ($endPeak+1) . "\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
			}
			elsif ($seq[$i] eq "C") {
				$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
			}
			elsif ($seq[$i] eq "G" and $i != 0 and defined $seq[$i-1] and $seq[$i-1] ne "C") {
				$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
				$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
			}
		}
		$head{CHH}{pos}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{pos}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t+\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
		$head{CHH}{neg}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CHH\"; transcript_id \"MYNUMBER.$gene.CHH\"\n";
		$head{CPG}{neg}{$gene} .= "$chr\tNA\texon\t$end0\t$end1\t0\t-\t0\tgene_id \"MYNUMBER.$gene.CPG\"; transcript_id \"MYNUMBER.$gene.CPG\"\n";
	}
#	foreach my $type (sort keys %head) {
#		next if $type eq "BLOCK";
#		foreach my $strand (sort keys %{$head{$type}}) {
#			print $GTFPosfile "track name=\"$type\_header +\" itemRgb=On\n" if $strand eq "pos";
#			print $GTFPosfile "track name=\"$type\_header +\" itemRgb=On\n" if $strand eq "pos";
#			print $GTFNegfile "track name=\"$type\_header -\" itemRgb=On\n" if $strand eq "neg";
#			print $GTFNegfile "track name=\"$type\_header -\" itemRgb=On\n" if $strand eq "neg";
#			foreach my $gene (sort keys %{$head{$type}{$strand}}) {
#				print $GTFPosfile "$head{$type}{$strand}{$gene}\n" if $strand eq "pos";
#				print $GTFNegfile "$head{$type}{$strand}{$gene}\n" if $strand eq "neg";
#			}
#		}
#	}
	close $faIn;
	return(\%head);
}


# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

sub parse {
	my ($gene, $file, $coor, $out) = @_;
	print "Parsing gene=$gene, file=$file\n";
	my $chr = $coor->{chr};
	my $beg = $coor->{beg};my ($beg0, $beg1) = ($beg, $beg + 1);
	my $end = $coor->{end};my ($end0, $end1) = ($end - 1, $end);
	my $strand = $coor->{strand};
	my $peak = 0;
	#hclust_R($file);# produces $file.Clust
	my %peak;
	open (my $in, "<", "$file") or die "Failed to read from $file.Clust: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		$line =~ s/["']//g;$line =~ s/^SEQ_//;
		my ($name, @vals) = split ("\t", $line);
		$name = "MYNUMBER.$gene.$name";
		my $begGTF = "$chr\tNA\texon\t$beg0\t$beg1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		my $endGTF = "$chr\tNA\texon\t$end0\t$end1\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"";
		$peak{$peak}{print} .= "$begGTF\n";
		my ($peakMinI, $peakMaxI, $peakMin, $peakMax) = (int($end/200+0.5), int($beg/200+0.5), $end, $beg);#to determine cluster
		for (my $i = 0;$i < @vals;$i++) {
			next if $vals[$i] =~ /^[62]$/;
			my $begPeak = $i + $beg + $bufferL;# 0 based
			my $endPeak = $begPeak + 1;# 1 based
			$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[9]$/;
#			$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[1]$/;
			$peak{$peak}{print} .= "$chr\tNA\tCDS\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[5]$/;
			#$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[4]$/;
			#$peak{$peak}{print} .= "$chr\tNA\texon\t$begPeak\t$endPeak\t0\t$strand\t0\tgene_id \"$name\"; transcript_id \"$name\"\n" if $vals[$i] =~ /^[0]$/;
			
			# For sorting
			my $begPeakIndex = int($begPeak/200+0.5);
			my $endPeakIndex = int($endPeak/200+0.5);
			if ($peakMinI > $begPeakIndex and $vals[$i] =~ /^[59]$/) {
				$peakMinI = $begPeakIndex;
				$peakMin  = $begPeak;
			}
			if ($peakMaxI < $endPeakIndex and $vals[$i] =~ /^[59]$/) {
				$peakMaxI = $endPeakIndex;
				$peakMax  = $endPeak;
			}
		}
		$peak{$peak}{print} .= "$endGTF\n";
		$peak{$peak}{begI} = $peakMinI;
		$peak{$peak}{endI} = $peakMaxI;
		$peak{$peak}{beg} = $peakMin;
		$peak{$peak}{end} = $peakMax;
		$peak++;
	}
	my $countz = 0;
	my $strand2 = $strand eq "+" ? "pos" : $strand eq "-" ? "neg" : die "Cannot determine strand ($strand), is not + or - !\n";
	foreach my $peak (sort {$peak{$b}{begI} <=> $peak{$a}{begI} || $peak{$b}{endI} <=> $peak{$a}{endI} || $peak{$b}{beg} <=> $peak{$a}{beg} || $peak{$b}{end} <=> $peak{$a}{end}} keys %peak) {
		my $countz2 = int($countz/100);
		my $PEAK = $peak{$peak}{print}; $PEAK =~ s/MYNUMBER/$countz2/g;
		print $out $PEAK;
		if ($countz % 100 == 0) {
			die "type = CHH, GENE = $gene, STRAND=$strand2\n" if not defined $head{CHH}{$strand2}{$gene};
			my $CHH = $head{CHH}{$strand2}{$gene}; $CHH =~ s/MYNUMBER/$countz2/g;
			my $CPG = $head{CPG}{$strand2}{$gene}; $CPG =~ s/MYNUMBER/$countz2/g;
			my $BLOCK = $head{BLOCK}{$strand2}{$gene}; $BLOCK =~ s/MYNUMBER/$countz2/g;
			print $out $CHH;
			print $out $CPG;
			print $out $BLOCK;
		}
		$countz ++;
	}
}

sub hclust_R {
	my ($reads) = @_;
	my $output = "$reads.Clust";
	open (my $out, ">", "$output.R");
	print $out "
	library(\"GMD\")
	df = read.table(\"./$reads\", sep=\"\t\", row.names=1)
	df2 = df
	df2[df2 != 9] = 0
	df\$sum = apply(df2,1,sum)
	df = df[order(-df\$sum),]
	dimz = dim(df)[1]
	df = df[1:dimz,]
	df = subset(df,select=-sum)
	df2 = df
	df2[df2 !=9] = 0
	h = heatmap.3(
		df2,
		dendrogram=\"row\",
		Rowv=TRUE, Colv=FALSE,
		labRow=FALSE, labCol=FALSE
	)
	id = h\$rowInd
	id = rev(id)
	df = df[id,]
	dev.off()
	write.table(df, file = \"./$output\", sep = \"\t\")
	";
	close $out;
	system("R --vanilla --no-save < $output.R") == 0 or die "Failed to run Rscript hclust ($output.R): $!\n";
	system("rm $output.R") == 0 or die "Failed to remove $output.R: $!\n";
}
