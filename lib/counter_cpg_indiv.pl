#!/usr/bin/perl

use strict; use warnings; use Cwd qw(abs_path); use File::Basename qw(dirname); use Getopt::Std;
use vars qw($opt_w $opt_s $opt_c $opt_m $opt_W $opt_f $opt_B $opt_a $opt_b $opt_M $opt_A $opt_o $opt_1 $opt_v);
getopts("o:Ww:s:cm:fBa:b:MA1v");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib;
use FAlite;
use Thread::Queue;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0);
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

my $date = getDate();

################
# ARGV Parsing #
###############

my ($input) = @ARGV;

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N [options] $LCY<fasta>$N

Options:
-1: include cpg dens and non-weighted atskew and gcskew
-o: outdir
-A: use all nucleotide instead of window length
-a: Nucleotide 1 (C)
-b: Nucletodie 2 (G)
    gcskew will be nucletode 2 vs 1 (GC gcskew)
-w: Window size [default: 200]
-s: Step size [default: 10]

";

die $usage unless @ARGV == 1;

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

my %outbuff;
my @outputtypes = qw(gccont gcwskew purineskew atwskew);
@outputtypes = qw(cpgdens gccont gcskew gcwskew purineskew atskew atwskew) if defined $opt_1;

foreach my $outputtype (@outputtypes) {
	my $outputfile = defined $opt_W ? "$outdir/$outName.$outputtype.wig" : "$outdir/$outName.$outputtype.tsv";
	open ($outbuff{$outputtype}, ">", $outputfile) or die "Cannot write to file=$outputfile type=$outputtype: $!\n";
}

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
	my ($res) = dinuc_window_NO_N($seq,$N1,$N2,$window,$step,$all);
	next if not defined $res or (defined $res and $res =~ /NULL/);
	foreach my $type (sort keys %{$res}) {
		if (not defined $opt_M) {
			$result{name}{$def}{$type} = $res->{$type};
			next;
		}
		for (my $i = 0; $i < @{$res->{$type}}; $i++) {
			push(@{$result{type}{$type}{$i}}, $res->{$type}[$i]);
		}
	}
}
my $total_read = (keys %result);
if (defined $opt_M) {
	print "DEFINED M\n";
	my $namez = $outName;
	($namez) = $outName =~ /^(.+)\.fa$/ if $outName =~ /\.fa$/;

	foreach my $type (sort keys %{$result{type}}) {
		foreach my $i (sort {$a <=> $b} keys %{$result{type}{$type}}) {
			my $totalz  = scalar(@{$result{type}{$type}{$i}});
			my $mean 	= mean(@{$result{type}{$type}{$i}});
			my $tmm 		= tmm(@{$result{type}{$type}{$i}});
			my $median 	= median(@{$result{type}{$type}{$i}});
			my $se 		= se(@{$result{type}{$type}{$i}});
			my $tmmse 	= tmmse(@{$result{type}{$type}{$i}});
			printf {$outbuff{$type}} "$i\t$totalz\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $mean,$tmm,$median,$se,$tmmse;
		}
	}
}
else {
	foreach my $name (sort keys %{$result{name}}) {
		my $start = int(0.5 + ($window / 2));
		foreach my $type (sort keys %{$result{name}{$name}}) {
			print {$outbuff{$type}} "fixedStep chrom=$name start=$start span=$step step=$step\n" if defined $opt_W;
			print {$outbuff{$type}} "$name" if not defined $opt_W;
		}
		
		foreach my $type (sort keys %{$result{name}{$name}}) {
			my $val = 0;
			if (defined($opt_m)) {
				$val = mean(@{$result{name}{$name}{$type}}) if $opt_m eq "mean" or $opt_m eq "ave";
				$val = median(@{$result{name}{$name}{$type}}) if $opt_m eq "med" or $opt_m eq "median";
				$val = min(@{$result{name}{$name}{$type}}) if $opt_m eq "min";
				$val = max(@{$result{name}{$name}{$type}}) if $opt_m eq "max";
				printf {$outbuff{$type}} "%.3f\n", $val;
			}
			else {
				for (my $i = 0; $i < @{$result{name}{$name}{$type}}; $i++) {
					$val = $result{name}{$name}{$type}[$i];
					printf {$outbuff{$type}} "\t%.3f", $val;
				}
				print {$outbuff{$type}} "\n";
			}
		}
	}
}
print "\nOutput:\n";
foreach my $outputtype (@outputtypes) {
	my $outputfile = defined $opt_W ? "$outdir/$outName.$outputtype.wig" : "$outdir/$outName.$outputtype.tsv";
	my ($md5) = `$md5script $outputfile` =~ /^(\w+) /;
	print "- $LGN$outputfile$N (md5=$YW$md5$N)\n";
}
print "\n";

my $name = getFilename($input);
if (defined $opt_B) {
	system("samtools faidx $input") == 0 or die "Failed to do samtools faidx $input: $!\n" if not -e "$input.fai";
	system("wigToBigWig -clip $outName.cpgdens.wig $input.fai $outName.cpgdens.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.cpgdens.wig $input.fai $outName.cpgdens.bw: $!\n";
	system("wigToBigWig -clip $outName.gccont.wig $input.fai $outName.gccont.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.gccont.wig $input.fai $outName.gccont.bw: $!\n";
	system("wigToBigWig -clip $outName.gcskew.wig $input.fai $outName.gcskew.bw") == 0 or die "Failed to do wigToBigWig -clip $outName.gcskew.wig $input.fai $outName.gcskew.bw: $!\n";
	system("aws s3 cp $outName.cpgdens.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.cpgdens.bw s3://muhucsc/: $!\n";
	system("aws s3 cp $outName.gccont.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.gccont.bw s3://muhucsc/: $!\n";
	system("aws s3 cp $outName.gcskew.bw s3://muhucsc/") == 0 or die "Failed to do aws s3 cp $outName.gcskew.bw s3://muhucsc/: $!\n";


	print "
track type=bigWig name=\"$outName CpG cpgdensity\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.cpgdens.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=0.2:1.2 windowingFunction=\"mean\" autoScale=off color=0,50,200
track type=bigWig name=\"$outName GC gccontent\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.gccont.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=0.3:0.8 windowingFunction=\"mean\" autoScale=off color=0,125,0
track type=bigWig name=\"$outName GC gcskew\" bigDataUrl=https://s3-us-west-1.amazonaws.com/muhucsc/$outName.gcskew.bw visibility=2 maxHeightPixels=40:40:40 viewLimits=-0.25:0.25 windowingFunction=\"mean\" autoScale=off color=200,50,0
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

	my $res;
   #my (@cpgdens, @gccont, @gcskew, @gcwskew, @purineskew, @atskew, @atwskew);
   my $count = -1;
   for (my $i = 0; $i <= length($sequence)-$window_size; $i += $step_size) { # change i to 1 for 1 bp window step
      $count++;

      my $seq_part = substr($sequence, $i, $window_size);
      my ($cpgdens_count, $nuc1_count, $nuc2_count) = (0, 0, 0);
		my ($Acount, $Ccount, $Gcount, $Tcount) = (0,0,0,0);
      # Count by Transliterate #
		$Acount = $seq_part =~ tr/A/A/;
		$Ccount = $seq_part =~ tr/C/C/;
		$Gcount = $seq_part =~ tr/G/G/;
		$Tcount = $seq_part =~ tr/T/T/;

      $nuc1_count = $nuc1 eq "A" ? $Acount : $nuc1 eq "T" ? $Tcount : $nuc1 eq "G" ? $Gcount : $Ccount;
      $nuc2_count = $nuc2 eq "A" ? $Acount : $nuc2 eq "T" ? $Tcount : $nuc2 eq "G" ? $Gcount : $Ccount;

      while ($seq_part =~ /$nuc1$nuc2/ig) {
         $cpgdens_count++;
      }
		if (defined $opt_1) {
	      $res->{cpgdens}[$count] = $nuc1_count * $nuc2_count 			== 0 ? 0 : ($cpgdens_count * $window_size) / ($nuc1_count * $nuc2_count);	
			$res->{atskew}[$count] 	= ($Acount+$Tcount) 			 			== 0 ? 0 : ($Acount-$Tcount)/($Acount+$Tcount);
	      $res->{gcskew}[$count] 	= $nuc1_count + $nuc2_count 			== 0 ? 0 : ($nuc2_count - $nuc1_count) / ($nuc2_count + $nuc1_count);
		}
      $res->{gccont}[$count] 		= ($nuc2_count + $nuc1_count) / $window_size;
      $res->{gcwskew}[$count] 	= $nuc1_count + $nuc2_count 			== 0 ? 0 : ($nuc2_count - $nuc1_count) / $window_size;
		$res->{purineskew}[$count] = ($Acount+$Gcount+$Ccount+$Tcount) == 0 ? 0 : ($Acount+$Gcount-$Ccount-$Tcount) / $window_size;
		$res->{atwskew}[$count] 	= ($Acount+$Tcount) 						== 0 ? 0 : ($Acount-$Tcount) / $window_size;
		
 
		foreach my $type (sort keys %{$res}) {
			die "fatal at sub DINUC_WINDOW_NO_N: type=$type not defined\n" if not defined($res->{$type}[$count]);
		}
   }
   return($res);
}
