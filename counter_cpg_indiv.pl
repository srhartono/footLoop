#!/usr/bin/perl

use strict; use warnings; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname); use Getopt::Std;
use vars qw($opt_w $opt_s $opt_c $opt_m $opt_W $opt_f $opt_B $opt_a $opt_b $opt_M $opt_A $opt_o);
getopts("o:Ww:s:cm:fBa:b:MA");


#########
# BEGIN #
#########

BEGIN {
   my ($bedtools) = `bedtools --version`;
   my ($bowtie2) = `bowtie2 --version`;
   my ($bismark) = `bismark --version`;
   my ($bismark_genome_preparation) = `bismark_genome_preparation --version`;
   print "\n\n\e[1;33m ------------- BEGIN ------------ \e[0m\n";
   if (not defined $bedtools or $bedtools =~ /command not found/ or $bedtools =~ /bedtools v?([01].\d+|2\.0[0-9]|2\.1[0-6])/) {
      print "Please install bedtools at least version 2.17 before proceeding!\n";
      $bedtools = 0;
   }
   print "\n- \e[1;32m bedtools v2.17+ exists:\e[0m " . `which bedtools` if $bedtools ne 0;
   die if $bedtools eq 0;
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
#  print "\n\n\e[1;33m ------------ BEGIN ------------> \e[0m\n";
}
use myFootLib; use FAlite;
use Thread::Queue;

################
# ARGV Parsing #
###############


my $date = getDate();


my ($input) = @ARGV;

die "\nUsage: $YW$0$N [options] $LCY<fasta>$N

Options:
-o: outdir
-A: count all nucleotide in the sequence
-a: Nucleotide 1 (C)
-b: Nucletodie 2 (G)
	Skew will be nucletode 2 vs 1 (GC skew)
-w: Window size [default: 200]
-s: Step size [default: 10]

" unless @ARGV == 1;

my ($outFolder, $outName) = getFilename($input, "folderfull") if not ($opt_c);
my $outdir = defined $opt_o ? $opt_o : $outFolder;
($outName) = $input =~ /\w+\_\w+\.(\w+)\.fa/ if ($opt_c);
die "Undefined name (maybe not a CEGMA format?\n" if ($opt_c) and not defined($outName);
my $N1 = defined $opt_a ? $opt_a : "C";
my $N2 = defined $opt_b ? $opt_b : "G";
my $window = defined($opt_w) ? $opt_w : 200;
my $step   = defined($opt_s) ? $opt_s : 10;
my $ext = defined $opt_W ? "wig" : "tsv";
open (my $in, "<", $input) or die;
my $outputDens = defined $opt_W ? "$outdir/$outName\.dens.wig" : "$outdir/$outName\.dens\.tsv";
my $outputCont = defined $opt_W ? "$outdir/$outName\.cont.wig" : "$outdir/$outName\.cont\.tsv";
my $outputSkew = defined $opt_W ? "$outdir/$outName\.skew.wig" : "$outdir/$outName\.skew\.tsv";
open (my $outDens, ">", $outputDens) or die "Cannot write to $outputDens: $!\n";
open (my $outCont, ">", $outputCont) or die "Cannot write to $outputCont: $!\n";
open (my $outSkew, ">", $outputSkew) or die "Cannot write to $outputSkew: $!\n";
my $Q = new Thread::Queue;
my $def;
my $record_number = 0;
my $fasta = new FAlite($in);
while (my $entry = $fasta->nextEntry) {
	my $def = $entry->def; $def =~ s/>//;
	if ($opt_f) {$def =~ s/ dna:.+$//; $def =~ s/^>//;}
	my $seq = $entry->seq;
	my @all = ($def, $seq, $record_number);
	$Q->enqueue(\@all);
	$record_number++;
}
$Q->end;
my $all = 1 if defined $opt_A;
my $count = 0;
my %result;
while ($Q->pending) {
	$count = $Q->pending();
	my ($def, $seq, $record_number) = @{$Q->dequeue};
	$seq = uc($seq);
	my ($dens, $cont, $skew) = dinuc_window_NO_N($seq,$N1,$N2,$window,$step,$all);
	next if $dens =~ /NULL/;
	$result{$record_number}{name} = $def;
	$result{$record_number}{dens}  = $dens;
	$result{$record_number}{cont}   = $cont;
	$result{$record_number}{skew} = $skew;
}

if (defined $opt_M) {
	my (@densSE, @contSE, @skewSE);
	my (@dens, @cont, @skew);
	my $total_read = (keys %result);
	foreach my $count (sort {$a <=> $b} keys %result) {
		my @dens1 = @{$result{$count}{dens}};
		my @cont1 = @{$result{$count}{cont}};
		my @skew1 = @{$result{$count}{skew}};
		# get mean
		for (my $i = 0; $i < @dens1; $i++) {
			$dens[$i] += $dens1[$i] eq "NA" ? 0 : ($dens1[$i] / $total_read);
			$cont[$i] += $cont1[$i] eq "NA" ? 0 : ($cont1[$i] / $total_read);
			$skew[$i] += $skew1[$i] eq "NA" ? 0 : ($skew1[$i] / $total_read);
		}
	}
	foreach my $count (sort {$a <=> $b} keys %result) {
		my @dens1 = @{$result{$count}{dens}};
		my @cont1 = @{$result{$count}{cont}};
		my @skew1 = @{$result{$count}{skew}};
		# get SD^2
		for (my $i = 0; $i < @dens1; $i++) {
			$densSE[$i] += $dens1[$i] eq "NA" ? 0 : (($dens1[$i] - $dens[$i])**2 / ($total_read-1));
			$contSE[$i] += $cont1[$i] eq "NA" ? 0 : (($cont1[$i] - $cont[$i])**2 / ($total_read-1));
			$skewSE[$i] += $skew1[$i] eq "NA" ? 0 : (($skew1[$i] - $skew[$i])**2 / ($total_read-1));
		}
	}
	my $namez = $outName;
	($namez) = $outName =~ /^(.+)\.fa$/ if $outName =~ /\.fa$/;
	# get SE  = sqrt(SD^2) / sqrt(total_read)
	for (my $i = 0; $i < @densSE; $i++) {
		my $densSD = sqrt($densSE[$i]);
		my $contSD = sqrt($contSE[$i]);
		my $skewSD = sqrt($skewSE[$i]);
		my $densSE = $densSD / sqrt($total_read);
		my $contSE = $contSD / sqrt($total_read);
		my $skewSE = $skewSD / sqrt($total_read);
		print $outDens "$i\t$dens[$i]\t$densSD\t$densSE\t$namez\n";
		print $outCont "$i\t$cont[$i]\t$contSD\t$contSE\t$namez\n";
		print $outSkew "$i\t$skew[$i]\t$skewSD\t$skewSE\t$namez\n";
	}
}
else {
foreach my $count (sort {$a <=> $b} keys %result) {
	my @dens1  = @{$result{$count}{dens}};
	my @cont1   = @{$result{$count}{cont}};
	my @skew1 = @{$result{$count}{skew}};
	my $name  = $result{$count}{name};
	my $start = int(0.5 + ($window / 2));
	if ($opt_W) {
		print $outDens "fixedStep chrom=$name start=$start span=$step step=$step\n";
		print $outCont "fixedStep chrom=$name start=$start span=$step step=$step\n";
		print $outSkew "fixedStep chrom=$name start=$start span=$step step=$step\n";
	}
	else {
		print $outDens "$name\t";
		print $outCont "$name\t";
		print $outSkew "$name\t";
		print "$YW$count$N: $name\t";
	}
	
	if (defined($opt_m)) {
		my ($dens, $cont, $skew) = (0,0,0);
		if ($opt_m eq "ave") {
			for (my $i = 0; $i < @dens1; $i ++) {
				$dens = $dens1[$i] eq "NA" ? $dens : $dens + ($dens1[$i] / @dens1);
				$cont = $cont1[$i] eq "NA" ? $cont : $cont + ($cont1[$i] / @cont1);
				$skew = $skew1[$i] eq "NA" ? $skew : $skew + ($skew1[$i] / @skew1);
			}
		}
		elsif ($opt_m eq "med") {
			my $count = 0;
			for (my $i = int(@dens1/3); $i <= int(@dens1*2/3); $i ++) {
				$count ++;
				$dens = $dens1[$i] eq "NA" ? $dens : $dens + $dens1[$i];
				$cont = $cont1[$i] eq "NA" ? $cont : $cont + $cont1[$i];
				$skew = $skew1[$i] eq "NA" ? $skew : $skew + $skew1[$i];
			}
			$dens /= $count;
			$cont /= $count;
			$skew /= $count;
		}
		elsif ($opt_m eq "max") {
			for (my $i = 0; $i < @dens1; $i ++) {
				$dens = $dens1[$i] if $dens1[$i] ne "NA" and $dens < $dens1[$i];
				$cont = $cont1[$i] if $cont1[$i] ne "NA" and $cont < $cont1[$i];
				$skew = $skew1[$i] if $skew1[$i] ne "NA" and $skew < $skew1[$i];
			}
		}

		elsif ($opt_m eq "max5skew") {
			my @skewsort;
			for (my $i = 0; $i < @skew1; $i++) {
				next if $skew1[$i] eq "NA";
				push(@skewsort, $skew1[$i]);
			}
			@skewsort = sort {$b <=> $a} (@skewsort);
			my $counts = 0;
			for (my $i = 0; $i < @skewsort; $i ++) {
				$dens = "NA" if $dens1[$i] ne "NA";
				$cont = "NA" if $cont1[$i] ne "NA";
				$skew += $skewsort[$i] if $skewsort[$i] ne "NA";
				$count ++;
				last if $i > int(0.05*@skewsort);
			}
			$skew /= $count;
		}
		printf $outDens "%.3f\n", $dens;
		printf $outCont "%.3f\n", $cont;
		printf $outSkew "%.3f\n", $skew;
	}
	else {
		for (my $i = 0; $i < @dens1; $i ++) {
			my $lineEnd = defined $opt_W ? "\n" : $i != @dens1 - 1 ? "\t" : "\n";
			my $densz = $dens1[$i] eq "NA" ? 0 : int($dens1[$i]*1000+0.5)/1000;
			my $contz = $cont1[$i] eq "NA" ? 0 : int($cont1[$i]*1000+0.5)/1000;
			my $skewz = $skew1[$i] eq "NA" ? 0 : int($skew1[$i]*1000+0.5)/1000;
			printf $outDens "$densz$lineEnd";
			printf $outCont "$contz$lineEnd";
			printf $outSkew "$skewz$lineEnd";
			print "$YW$count$N:$LGN$i$N: $skewz, lineEnd=$LGN$lineEnd$N\n";
		}
	}
}
}
print "\nOutput:
- $LCY$outputDens$N
- $LGN$outputCont$N
- $LRD$outputSkew$N
\n";

my $name = getFilename($input);
if (defined $opt_B) {
	system("samtools faidx $input") == 0 or die "Failed to do samtools faidx $input: $!\n" if not -e "$input.fai";
	system("wigToBigWig -clip $outName.dens.wig $input.fai $outName.dens.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.dens.wig $input.fai $outName.dens.bw: $!\n";
	system("wigToBigWig -clip $outName.cont.wig $input.fai $outName.cont.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.cont.wig $input.fai $outName.cont.bw: $!\n";
	system("wigToBigWig -clip $outName.skew.wig $input.fai $outName.skew.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.skew.wig $input.fai $outName.skew.bw: $!\n";
	system("aws s3 cp $outName.dens.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.dens.bw s3://muhucsc/: $!\n";
	system("aws s3 cp $outName.cont.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.cont.bw s3://muhucsc/: $!\n";
	system("aws s3 cp $outName.skew.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.skew.bw s3://muhucsc/: $!\n";


	print "
track type=bigWig name=\"$outName CpG Density\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.dens.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=0.2:1.2 windowingFunction=\"mean\" autoScale=off color=0,50,200
track type=bigWig name=\"$outName GC Content\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.cont.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=0.3:0.8 windowingFunction=\"mean\" autoScale=off color=0,125,0
track type=bigWig name=\"$outName GC Skew\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.skew.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=-0.25:0.25 windowingFunction=\"mean\" autoScale=off color=200,50,0
";
}


sub dinuc_window_NO_N {
   my ($sequence, $nuc1, $nuc2, $window_size, $step_size, $all) = @_;
   return "NULL" if not defined($sequence);
   $window_size = 100 if not defined($window_size);
   $step_size = 1 if not defined($step_size);

   # RETURN NULL IF NUMBER OF N IS MORE THAN 10% OF GENE LENGTH
   my $N_Count = $sequence =~ tr/Nn/Nn/;
   my $length_Sequence = length($sequence);
   return "NULL" if $N_Count >= 0.1 * $length_Sequence and not defined $all;
   $window_size = length($sequence) if defined $all;


   my (@dens, @cont, @skew, @wskew);
   my $count = -1;
   for (my $i = 0; $i <= length($sequence)-$window_size; $i += $step_size) { # change i to 1 for 1 bp window step
      $count++;

      my $seq_part = substr($sequence, $i, $window_size);
      my ($dens_count, $nuc1_count, $nuc2_count) = (0, 0, 0);

      # Count by Transliterate #
      $nuc1_count = $seq_part =~ tr/A/A/ if $nuc1 =~ /A/;
      $nuc1_count = $seq_part =~ tr/T/T/ if $nuc1 =~ /T/;
      $nuc1_count = $seq_part =~ tr/G/G/ if $nuc1 =~ /G/;
      $nuc1_count = $seq_part =~ tr/C/C/ if $nuc1 =~ /C/;

      $nuc2_count = $seq_part =~ tr/A/A/ if $nuc2 =~ /A/;
      $nuc2_count = $seq_part =~ tr/T/T/ if $nuc2 =~ /T/;
      $nuc2_count = $seq_part =~ tr/G/G/ if $nuc2 =~ /G/;
      $nuc2_count = $seq_part =~ tr/C/C/ if $nuc2 =~ /C/;

      while ($seq_part =~ /$nuc1$nuc2/g) {
         $dens_count++;
      }
      my $dens = $nuc1_count * $nuc2_count == 0 ? 0 : ($dens_count * $window_size) / ($nuc1_count * $nuc2_count);
      my $cont = ($nuc2_count + $nuc1_count) / $window_size;
      my $skew = $nuc1_count + $nuc2_count == 0 ? 0 : ($nuc2_count - $nuc1_count) / ($nuc2_count + $nuc1_count);
      my $equalizer = $nuc2_count > $nuc1_count ? $nuc2_count : $nuc1_count;
      my $wskew = $nuc1_count + $nuc2_count == 0 ? 0 : (($nuc2_count - $nuc1_count) * ($equalizer / $window_size)) / ($nuc2_count + $nuc1_count);
      $dens[$count] = $seq_part =~ /^[Nn]+$/ ? "NA" : $dens;
      $cont[$count] = $seq_part =~ /^[Nn]+$/ ? "NA" : $cont;
      $skew[$count] = $seq_part =~ /^[Nn]+$/ ? "NA" : $skew;
      $wskew[$count] = $seq_part =~ /^[Nn]+$/ ? "NA" : $wskew;
      die "fatal at sub dinuc_window_NO_N: cpg dens not defined\n" if not defined($dens);
      die "fatal at sub dinuc_window_NO_N: gc content not defined\n" if not defined($cont);
      die "fatal at sub dinuc_window_NO_N: gc skew not defined\n" if not defined($skew);
      die "fatal at sub dinuc_window_NO_N: weighted gc skew not defined\n" if not defined($wskew);

   }
   return(\@dens, \@cont, \@skew, \@wskew);
}
