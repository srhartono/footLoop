#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use vars qw($opt_v $opt_i $opt_n $opt_m $opt_r $opt_s $opt_d);
getopts("vi:n:m:r:s:d:");

$opt_d = 98 if not defined $opt_d;
#$opt_s = 4200 if not defined $opt_s;
$opt_r = 20 if not defined $opt_r;
$opt_n = 250 if not defined $opt_n;
$opt_m = "2p" if not defined $opt_m;

my ($input1, $simtot, $MAXEDITS, $random, $minidentity) = ($opt_i, $opt_n, $opt_m, $opt_r, $opt_d);
die "
Usage: $0 -i <fastq.fq>

-s : [positive integer, default 4200] random number generator
-r : [positive integer, default: 20] total number of reads to be randomly chosen is up to this.
-n : [positive integer, default: 250] pbsim total read depth (total simulated read ~ depth x 4)
     Don't use too few (e.g. 20) as it'll create too few samples to be statistically meaningful
     The more -n, the longer it'll take (exponentially).
-d : [default: 98] percent identity for dedupe.sh minidentity parameter.
-m : [default: 2\% read length]  dedupe.sh maxedits parameter.
     put \"p\" after number to make it based on current read's length,
     otherwise, it'll directly translate into basepair. Examples:
     -m 2p: 2\% of read length
     -m 4p: 4\% of read length
     -m 2: 2 basepairs
     -m 10: 10 basepairs
     

Ideally 100\% simulated reads should either be duplicates (100\% identical) or containment (similar enough to be duplicates based on maxedits and minidentity)
However, even with maximum settings, only 85-90\% maximum reads are marked as duplicates/containment.
Therefore the goal is to find parameters that'll call 85-90\% reads as duplicates/containments, 
leaving 10-15\% remaining reads (false negative is at most 10-15\%).

On containment parametes:
- minidentity specifies how many mismatches to still be considered \"similar enough\".
  Default is 98 \%, as expected pacbio ccs mismatches is 2\% at most.
- maxedits specifies how many indels to still be considered \"similar enough\".
  Default -m is 2\% of read length, as expected pacbio ccs indels is 2\% at most.
  2\% is also the threshold which consistently call maximum \% simulated reads into containments
  (leaving 10-15\% \"duplicates\"). 
  More than 2\% will cause no/minuscule change on \% containment in simulated reads, 
  but on the real data it'll significantly call more reads into containment (false positives)
  Less than that it'll sometimes cause significantly less \% containment in simulated reads, hence false negatives.
  1-2\% causes some reads to have only 60-70\% containments from testing.


  \e[1;32mTherefore, we should use maxedits 2\% (more false positives) to be on the more conservative side, making sure
  that our actual data only have at most 10-15\% duplicates (false negatives).\e[0m

" unless defined $input1;

die "-s seed for random number generator must be positive integer (current: $opt_s)\n" if defined $opt_s and $opt_s !~ /^\d+/;
die "-d minidentity must be positive integer (current: $opt_d)\n" if $opt_d !~ /^\d+/;
die "-n total simulation must be positive integer (current: $opt_n)\n" if $opt_n !~ /^\d+/;
die "-m max edits must be positive integer with or w/o p (current: $opt_m)\n" if $opt_m !~ /^\d+\.?\d*p?$/;
die "-r total randomly chosen read must be positive integer (current: $opt_r)\n" if $opt_r !~ /^\d+e?\d*$/;
#die "-r must be <= 1000! (current: $opt_r)\n" if $opt_r > 1000;
die "-d must be <= 100! (current: $opt_d)\n" if $opt_d > 100;
srand($opt_s) if defined $opt_s;
my $rand = defined $opt_s ? "-s $opt_s" : "";
my $maxeditscheck = $MAXEDITS . "bp";
if ($MAXEDITS =~ /^\d+\.?\d*p$/) {
	($maxeditscheck) = $MAXEDITS =~ /^(\d+\.?\d*)p$/;
	die "-m maxedits is with p ($MAXEDITS) but number is greater than 100... ($maxeditscheck)\n" if $maxeditscheck > 100;
	$maxeditscheck = "$maxeditscheck \% of each read length\n";
}

if ($opt_r =~ /^\d+e\d+$/) {
	my $zero = 0;
	($random, $zero) = $opt_r =~ /^(\d+)e*(\d+)$/;
	$random = $random * 10**$zero;
}

print "\n\e[1;33m================================== RUN PARAMETES ======================================\e[0m\n\n";
print "

Run script             = $0 -i $opt_i -r $opt_r -n $opt_n -d $opt_d -m $opt_m $rand
-i Fastq file          = $input1
-r random chosen read  = $random
-n total simulation    = $simtot
-m maxedits            = $maxeditscheck

";

my ($total_linecount) = `zcat $input1 | wc -l` =~ /^(\d+)/ if $input1 =~ /.gz$/;
($total_linecount) = `wc -l $input1` =~ /^(\d+)/ if $input1 !~ /.gz$/;
die if $total_linecount % 4 != 0;
$total_linecount /= 4;
my @random; my %random;
for (my $i = 0; $i < $total_linecount; $i++) {
	push(@random, $i);
}
for (my $i = 0; $i < $total_linecount * 5; $i++) {
	my $rand1 = int(rand(@random));
	my $randval1 = $random[$rand1];
	my $rand2 = int(rand(@random));
	my $randval2 = $random[$rand1];
	$random[$rand1] = $randval2;
	$random[$rand2] = $randval1;
}
for (my $i = 0; $i < $total_linecount; $i++) {
	$random{$random[$i]} = 1;
	last if $i > $random;
}
print "Done randoming $random numbers!\n";

my %data; my $total_done = 0; my $linecount = -1; my $too_short = 0;
my $in1;
system("mkdir \"$input1.pbsim\"") == 0 or die "Failed to create directory $input1.pbsim: $!\n" if not -d "$input1.pbsim";
open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /gz$/;
open ($in1, "zcat $input1|") or die "Cannot read from zcat $input1: $!\n" if $input1 =~ /gz$/;
open (my $outStat, ">", "$input1.dedupe_stat.csv") or die;
print $outStat "no,read_number,length_seq,maxedits_bp,total_reads,perc_duplicates,perc_containment,perc_remaining\n";
while (my $line = <$in1>) {
	$linecount ++;
	chomp($line);
	my ($read) = $line; $read =~ s/^\@//;
	my $read2 = $read; 
	if ($read =~ /\/\d+\/ccs/) {
		($read2) = $read =~ /\/(\d+)\/ccs/;
	}
	elsif ($read =~ /^(\d+)\..+/) {
		($read2) = $read =~ /^(\d+)\./;
	}
	elsif ($read =~ /^.+\.(\d+)$/) {
		($read2) = $read =~ /^.+\.(\d+)$/;
	}
	else {
		$read2 = $read;
		$read2 =~ s/[\[\]\/\\\-\_\+\.\,\;\:\'\"\~\`\+\=\)\(\*\&\^\%\$\#\@\!]+//g;
	}
	die "read $read defined twice..\n" if defined $data{$read2};
	$line = <$in1>; chomp($line);
	my ($seq) = $line;
	my $lenseq = length($seq);
	my $maxedits = $MAXEDITS;
	if ($MAXEDITS =~ /^\d+\.?\d*p$/) {
		($maxedits) = $MAXEDITS =~ /^(\d+\.?\d*)p$/;
		$maxedits = int($maxedits/100*$lenseq+0.5);#100;#50;#int(0.02*$lenseq+0.5);
	}
	$line = <$in1>; chomp($line);
	#+
	$line = <$in1>; chomp($line);
	#qual
	next if not defined $random{$linecount};
	if ($lenseq < 1000) {
		$too_short ++;
#		print "\tlinecount=$linecount read number=$read2 too short ($lenseq), replaced with: ";
		my $add = 0;
		while ($linecount+$add < $total_linecount) {
			$add ++;
			next if defined $random{$linecount+$add};
			$random{$linecount+$add} = 2; $too_short --; 
			#print "linecount=" . ($linecount+$add). "\n"; 
			last;
		}
		if ($linecount + $add >= $total_linecount) {
			#print " ... Can't find replacement!\n";
		}
		next;
	}

	system("mkdir \"$input1.pbsim/$read2\"") == 0 or die "Failed to create directory $input1.pbsim/$read2: $!\n" if not -d "$input1.pbsim/$read2";
	open (my $out1, ">", "$input1.pbsim/$read2/$read2.fa") or die "Cannot write to $input1.pbsim/$read2/$read2.fa: $!\n";
	print $out1 ">$read2\n$seq\n";
	close $out1;
	

	print "\n\e[1;33m==================================== EXAMPLE ========================================\e[0m\n\n\n" if $total_done == 0;
	print "\e[1;33mEXAMPLE PIPELINE RAN ON \e[1;35m$read\e[0m (read number \e[1;35m$read2\e[0m)\n\n" if $total_done == 0;
	#run pbsim
	my $PBSIMCMD = "cd $input1.pbsim/$read2/ && pbsim --data-type CCS --model_qc /home/mitochi/src/PBSIM-PacBio-Simulator/data/model_qc_ccs --seed 1 --depth $simtot $read2.fa > /dev/null 2>&1";
	print "\t1. Example PBSIM: $PBSIMCMD\n" if $total_done == 0;
	system("$PBSIMCMD") == 0 or die "Failed to run $PBSIMCMD: $!\n";

	#run fq2fa
	my $FQ2FACMD = "fq2fa.pl $input1.pbsim/$read2/sd_0001.fastq > /dev/null 2>&1";
	print "\t2. Example fastq2fa: $FQ2FACMD\n" if $total_done == 0;
	system("$FQ2FACMD") == 0 or die "Failed to run $FQ2FACMD: $!\n";

	#run dedupe
	my $DEDUPECMD = "dedupe.sh ow minoverlap=5000 minidentity=$minidentity maxedits=$maxedits in=$input1.pbsim/$read2/sd_0001.fastq.fa out=$input1.pbsim/$read2/sd_0001.fastq.fa.rmdup > $input1.pbsim/$read2/RES.txt 2>&1";
	print "\t3. Example DEDUPE: $DEDUPECMD\n" if $total_done == 0;
	system("$DEDUPECMD") == 0 or die "Failed to run $DEDUPECMD: $!\n";

	#parse the result
	my @res = `cat $input1.pbsim/$read2/RES.txt`;
	print "\t4. Example resulting dedupe (\e[1;32m$input1.pbsim/$read2/RES.txt\e[0m)\n\n" . `cat $input1.pbsim/$read2/RES.txt` . "\n\n" if $total_done == 0;
	my ($tot, $dup, $con, $res) = (0,0,0,0);

	foreach my $line2 (@res) {
		chomp($line2);
		if ($line2 =~ /Input:\s+\d+ reads/) {
			($tot) = $line2 =~ /Input:\s+(\d+) reads/;
		}
		if ($line2 =~ /Duplicates:\s+\d+ reads/) {
			($dup) = $line2 =~ /Duplicates:\s+(\d+) reads/;
		}
		if ($line2 =~ /Containments:\s+\d+ reads/) {
			($con) = $line2 =~ /Containments:\s+(\d+) reads/;
		}
		if ($line2 =~ /Result:\s+\d+ reads/) {
			($res) = $line2 =~ /Result:\s+(\d+) reads/;
		}
	}
	$dup = $tot == 0 ? 0 : int(1000*$dup/$tot+0.5)/10;
	$con = $tot == 0 ? 0 : int(1000*$con/$tot+0.5)/10;
	$res = $tot == 0 ? 0 : int(1000*$res/$tot+0.5)/10;
	print "\n\e[1;33m====================================== RESULT =======================================\e[0m\n\n\n" if $total_done == 0;
	print "(also created in .csv format (google sheet friendly) at \e[1;33m$input1.dedupe_stat.csv\e[0m)\n\n" if $total_done == 0;
	print "\e[1;33m$total_done\e[0m. line number=$linecount ccs read_number=\e[1;33m$read2\e[0m; length_seq=\e[1;32m$lenseq\e[0m; maxedits=\e[1;32m$maxedits\e[0m; total_reads=\e[1;32m$tot\e[0m, duplicates=\e[1;32m$dup\%\e[0m, containment=\e[1;32m$con\%\e[0m, remaining=\e[1;32m$res\%\e[0m\n";
	print $outStat "$total_done,$read2,$lenseq,$maxedits,$tot,$dup,$con,$res\n";
	push(@{$data{tot}}, $tot);
	push(@{$data{dup}}, $dup);
	push(@{$data{con}}, $con);
	push(@{$data{res}}, $res);
	push(@{$data{lenseq}}, $lenseq);
	$total_done ++;
#	die "RES: tot=$tot dup=$dup con=$con res=$res\n" if
#	print "Done $total_done\n" if $total_done % 10 == 0;
	last if $total_done >= $random;
}
close $in1;
my $tot = int(100*tmm(@{$data{tot}})+0.5)/100;
my $totsd = int(100*tmmsd(@{$data{tot}})+0.5)/100;
my $dup = int(100*tmm(@{$data{dup}})+0.5)/100;
my $dupsd = int(100*tmmsd(@{$data{dup}})+0.5)/100;
my $con = int(100*tmm(@{$data{con}})+0.5)/100;
my $consd = int(100*tmmsd(@{$data{con}})+0.5)/100;
my $res = int(100*tmm(@{$data{res}})+0.5)/100;
my $ressd = int(100*tmmsd(@{$data{res}})+0.5)/100;
my $lenseq = int(100*tmm(@{$data{lenseq}})+0.5)/100;
my $lenseqsd = int(100*tmmsd(@{$data{lenseq}})+0.5)/100;

print "\n\e[1;33m================================== STATISTICS OF ALL ===================================\e[0m\n\n\n";
print "Out of \e[1;32m $total_done \e[0m randomly chosen reads (~ (4 x $simtot) simulated reads each):\n";
print "length seq  = \e[1;32m$lenseq\e[0m bp +/= \e[1;32m$lenseqsd \e[0m bp\n";
print "total reads = \e[1;32m$tot\e[0m reads +/= \e[1;32m$totsd\e[0m reads\n";
print "duplicates  = \e[1;32m$dup\e[0m \% +/= \e[1;32m$dupsd \e[0m \%\n";
print "containment = \e[1;32m$con\e[0m \% +/= \e[1;32m$consd \e[0m \%\n";
print "remaining   = \e[1;32m$res\e[0m \% +/= \e[1;32m$ressd \e[0m \%\n";
print "\n\e[1;33m========================================= DONE ====== ===================================\e[0m\n\n\n";

print "Output statistics in csv (google sheet friendly): \e[1;33m$input1.dedupe_stat.csv\e[0m\n\n";
sub tmm {
   my (@value) = @_;
   return(0) if @value == 0;
   my $mean = 0;
   @value = sort {$a <=> $b} @value;
   my $perc = int(0.05*@value);
   for (my $i = $perc; $i < @value-$perc; $i++) {
      $mean += $value[$i] / (@value-2*$perc);
   }
   return($mean);
}

sub mean {
       my (@value) = @_;
   return(0) if @value == 0;
   my $mean = 0;
   for (my $i =0 ; $i < @value; $i++) {
      $mean += $value[$i] / @value;
   }
        return($mean);
}
sub tmmsd {
        my (@value) = @_;
   my $mean = tmm(@value);
   return(0) if @value <= 1;
   my $total = @value;
   my $sd = 0;
   #print "VALUE @value\nMEAN $mean\n";
   @value = sort {$a <=> $b} @value;
   my $perc = int(0.05*@value+0.5);
   for (my $i = $perc; $i < @value-$perc; $i++) {
      $sd += (($value[$i] - $mean)**2) / (@value-2*$perc-1);
   }
   $sd = sqrt($sd);
   return($sd);
}
sub sd {
        my (@value) = @_;
   my $mean = mean(@value);
   return(0) if @value <= 1;
   my $total = @value;
   $total = 1 if $total == 0;
   my $sd = 0;
   #print "VALUE @value\nMEAN $mean\n";
   for (my $i = 0; $i < @value; $i++) {
      $sd += (($value[$i] - $mean)**2) / (@value-1);
   }
   $sd = sqrt($sd);
   return($sd);
}
