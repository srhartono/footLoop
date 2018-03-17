#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_n $opt_e $opt_r $opt_s $opt_m $opt_p);
getopts("vi:n:m:r:s:d:p:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   @INC = ($libPath, @INC);
}
use myFootLib;

$opt_m = 98 if not defined $opt_m;
#$opt_s = 4200 if not defined $opt_s;
$opt_r = 20 if not defined $opt_r;
$opt_n = 250 if not defined $opt_n;
$opt_e = "30" if not defined $opt_e;

my ($input1, $simtot, $random) = ($opt_i, $opt_n, $opt_r);


die "
Usage: $0 -i <fastq.fq>

-p : dedupe.sh additional parameter, COMMA separated. (E.g. mo=5000,k=30). 

     - \e[0;36mBy default, this script will run dedupe.sh with these parameters:
       \e[0;32m ow=t,k=31,nam=4,mo=200,e=30,mid=98\e[0m

     - Any parameter defined using -p will add/change the values of parameters above.
       Example, to use e=50, k=30, and c=t (maxedits 30, kmer 30, enabling cluster): 
       \e[0;32m-p e=50,k=30,c=t
       Results in: 
       \e[0;32m ow=t,\e[0;33mk=30\e[0;32m,nam=4,mo=200,\e[0;33me=50\e[0;32m,mid=98,\e[0;33mc=t\e[0m

     - Put 'NA' as value to not use it when running dedupe (thus using dedupe's default parameter value)
       E.g. to use e=50, k=30, c=t, but use dedupe's default mid (100): 
       \e[0;32m-p e=50,k=30,c=t,mid=NA\e[0m

     - on maxedits (e): 
     put \"p\" after number to make it based on current read's length,
     otherwise, it'll directly translate into basepair. Examples:
     e=2p: 2\% of read length
     e=4p: 4\% of read length
     e=2: 2 basepairs
     e=10: 10 basepairs

-s : [positive integer, default 4200] random number generator
-r : [positive integer, default: 20] total number of reads to be randomly chosen is up to this.
-n : [positive integer, default: 250] pbsim total read depth (total simulated read ~ depth x 4)
     Don't use too few (e.g. 20) as it'll create too few samples to be statistically meaningful
     The more -n, the longer it'll take (exponentially).

" unless defined $input1;

my $date = getDate();
my $logFile = "$input1.$date.logFile.txt";
open (my $outLog, ">", $logFile) or die "Failed to create footLoop_simtest.pl logFile: $!\n";

DIELOG($outLog, "-s seed for random number generator must be positive integer (current: $opt_s)\n") if defined $opt_s and $opt_s !~ /^\d+/;
DIELOG($outLog, "-n total simulation must be positive integer (current: $opt_n)\n") if $opt_n !~ /^\d+/;
#DIELOG($outLog, " max edits must be positive integer with or w/o p (current: $opt_e)\n" if $opt_e ne "NA") and $opt_e !~ /^\d+\.?\d*p?$/;
DIELOG($outLog, "-r total randomly chosen read must be positive integer (current: $opt_r)\n") if $opt_r !~ /^\d+e?\d*$/;
#DIELOG($outLog, "-r must be <= 1000! (current: $opt_r)\n") if $opt_r > 1000;
DIELOG($outLog, "-m must be <= 100! (current: $opt_m)\n") if $opt_m > 100;
srand($opt_s) if defined $opt_s;
my $rand = defined $opt_s ? "-s $opt_s" : "";

if ($opt_r =~ /^\d+e\d+$/) {
	my $zero = 0;
	($random, $zero) = $opt_r =~ /^(\d+)e*(\d+)$/;
	$random = $random * 10**$zero;
}

my ($ow, $k, $nam, $mo, $e, $mid) = ("t", 31, 4, 200, 30, 98);
my %param = ("ow"=>$ow,"k"=>$k,"nam"=>$nam,"mo"=>$mo,"e"=>$e,"mid"=>$mid);
if (defined $opt_p) {
	my @dedupe_param = split(",", $opt_p);
	foreach my $params (@dedupe_param) {
		my ($param, $value) = $params =~ /^(.+)=(.+)$/;
		DIELOG($outLog, "Undefined both param and value in params: \e[0;32m$params\e[0m\n") if not defined $param and not defined $value;
		DIELOG($outLog, "Undefined param (value=$value) in params: \e[0;32m$params\e[0m\n") if not defined $param;
		DIELOG($outLog, "Undefined value (param=$param) in params: \e[0;32m$params\e[0m\n") if not defined $value;
		$param = "mo" if $param eq "minoverlap";
		$param = "nam" if $param eq "numaffixmaps";
		$param = "e" if $param eq "maxedits";
		$param = "mid" if $param eq "minidentity";
		if ($value eq "NA") {
			delete $param{$param};
		}
		else {
			$param{$param} = $value;
		}
	}
}
my $dedupe_param;
foreach my $param (sort keys %param) {
	$dedupe_param .= " $param=$param{$param}";
}

my $MAXEDITS = defined $param{e} ? $param{e} : 0;
my $maxeditscheck = $MAXEDITS . "bp";
if ($MAXEDITS =~ /^\d+\.?\d*p$/) {
	($maxeditscheck) = $MAXEDITS =~ /^(\d+\.?\d*)p$/;
	DIELOG($outLog, "maxedits (e) is specified as percent ($MAXEDITS) but number is greater than 100... ($maxeditscheck)\n") if $maxeditscheck > 100;
	DIELOG($outLog, "maxedits (e) must be positive integer between 0-100 (current: $maxeditscheck)\n") if $maxeditscheck =~ /^\d+\.?\d*$/ and ($maxeditscheck < 0 or $maxeditscheck > 100);
	$maxeditscheck = "$maxeditscheck \% of each read length\n";
}
else {
	DIELOG($outLog, "maxedits (e) must be positive integer between 0-100 (current: $MAXEDITS)\n") if $MAXEDITS !~ /^\d+$/ or ($MAXEDITS =~ /^\d+$/ and ($MAXEDITS < 0 or $MAXEDITS > 100));
}
my $minidentity = defined $param{mid} ? $param{mid} : 100;
DIELOG($outLog, "minidentity (m) must be positive integer <= 100 (current: $minidentity)\n") if $minidentity !~ /^\d+$/ or ($minidentity =~ /^\d+$/ and $minidentity > 100);

#my $dedupe_param = defined $opt_p ? $opt_p : 
#my $dedupe_param = "ow=t,e=30,k=31,nam=4,minoverlap=200,maxedits=MAXEDITS,mid=MINIDENTITY";
#$dedupe_param =~ s/,/ /g;

my $dedupe_param_print = $dedupe_param;
$dedupe_param_print =~ s/ /\n/g;

LOG($outLog, "\n\e[1;33m================================== RUN PARAMETES ======================================\e[0m\n\n");
LOG($outLog, "

>This script parameters:
Run script             = $0 -i $opt_i -r $opt_r -n $opt_n $rand -p $dedupe_param
-i Fastq file          = $input1
-r random chosen read  = $random
-n total simulation    = $simtot

>dedupe.sh parameters:$dedupe_param_print

");

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
LOG($outLog, "Done randoming $random numbers!\n");

my %data; my $total_done = 0; my $linecount = -1; my $too_short = 0;
my $in1;
system("mkdir \"$input1.pbsim\"") == 0 or DIELOG($outLog, "Failed to create directory $input1.pbsim: $!\n") if not -d "$input1.pbsim";
$input1 =~ /.gz$/ ? LOG($outLog, "\nInput is gz\n") : LOG($outLog, "\nInput is NOT zipped file!\n");
open ($in1, "<", $input1) or DIELOG($outLog, "Cannot read from $input1: $!\n") if $input1 !~ /gz$/;
open ($in1, "zcat $input1|") or DIELOG($outLog, "Cannot read from zcat $input1: $!\n") if $input1 =~ /gz$/;
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
	DIELOG($outLog, "read $read defined twice..\n") if defined $data{$read2};
	$line = <$in1>; chomp($line);
	my ($seq) = $line;
	my $lenseq = length($seq);

	$line = <$in1>; chomp($line);
	#+

	$line = <$in1>; chomp($line);
	#qual

	next if not defined $random{$linecount};
	if ($lenseq < 1000) {
		$too_short ++;
#		LOG($outLog, " linecount=$linecount read number=$read2 too short ($lenseq), replaced with: ");
		my $add = 0;
		while ($linecount+$add < $total_linecount) {
			$add ++;
			next if defined $random{$linecount+$add};
			$random{$linecount+$add} = 2; $too_short --; 
			#LOG($outLog, "linecount=" . ($linecount+$add). "\n"); 
			last;
		}
		if ($linecount + $add >= $total_linecount) {
			#LOG($outLog, " ... Can't find replacement!\n");
		}
		next;
	}

	my $dedupe_param2 = $dedupe_param;
	my $maxedits = $MAXEDITS;
	if ($MAXEDITS =~ /^\d+\.?\d*p$/) {
		($maxedits) = $MAXEDITS =~ /^(\d+\.?\d*)p$/;
		LOG($outLog, "maxedits = $maxedits / 100 * $lenseq = " . ($maxedits/100*$lenseq) . "\n");
		$maxedits = int($maxedits/100*$lenseq+0.5);#100;#50;#int(0.02*$lenseq+0.5);
		$dedupe_param2 =~ s/e=\d+p?/e=$maxedits/;
	}

	system("mkdir \"$input1.pbsim/$read2\"") == 0 or DIELOG($outLog, "Failed to create directory $input1.pbsim/$read2: $!\n") if not -d "$input1.pbsim/$read2";
	open (my $out1, ">", "$input1.pbsim/$read2/$read2.fa") or DIELOG($outLog, "Cannot write to $input1.pbsim/$read2/$read2.fa: $!\n");
	print $out1 ">$read2\n$seq\n";
	close $out1;
	

	LOG($outLog, "\n\e[1;33m==================================== EXAMPLE ========================================\e[0m\n\n\n") if $total_done == 0;
	LOG($outLog, "\e[1;33mEXAMPLE PIPELINE RAN ON \e[1;35m$read\e[0m (read number \e[1;35m$read2\e[0m)\n\n") if $total_done == 0;
	#run pbsim
	my $PBSIMCMD = "cd $input1.pbsim/$read2/ && pbsim --data-type CCS --model_qc /home/mitochi/src/PBSIM-PacBio-Simulator/data/model_qc_ccs --seed 1 --depth $simtot $read2.fa > /dev/null 2>&1";
	LOG($outLog, "\t1. Example PBSIM\n\t\e[0;32m :: $PBSIMCMD :: \e[0m\n\n") if $total_done == 0;
	system("$PBSIMCMD") == 0 or DIELOG($outLog, "Failed to run $PBSIMCMD: $!\n");

	#run fq2fa
	my $FQ2FACMD = "cd $input1.pbsim/$read2/ && fq2fa.pl sd_0001.fastq > /dev/null 2>&1";
	LOG($outLog, "\t2. Example fq2fa:\n\t\e[0;32m :: $FQ2FACMD :: \e[0m\n\n") if $total_done == 0;
	system("$FQ2FACMD") == 0 or DIELOG($outLog, "Failed to run $FQ2FACMD: $!\n");

	#run dedupe
#	my $DEDUPECMD = "cd $input1.pbsim/$read2/ && dedupe.sh ow e=30 mid=98 k=31 nam=4 minoverlap=200 numaffixmaps=4 minidentity=$minidentity maxedits=$maxedits in=sd_0001.fastq.fa out=sd_0001.fastq.fa.rmdup > $read2.RES.txt 2>&1";
	my $DEDUPECMD = "cd $input1.pbsim/$read2/ && dedupe.sh $dedupe_param2 in=sd_0001.fastq.fa out=sd_0001.fastq.fa.rmdup > $read2.RES.txt 2>&1";
	LOG($outLog, "\t3. Example DEDUPE:\n\t\e[0;32m :: $DEDUPECMD :: \e[0m\n\n") if $total_done == 0;
	system("$DEDUPECMD") == 0 or DIELOG($outLog, "Failed to run $DEDUPECMD: $!\n");

	#parse the result
	my @res = `cat $input1.pbsim/$read2/$read2.RES.txt`;
	if ($total_done == 0) {
		LOG($outLog, "\t4. Example resulting dedupe (\e[1;32m$input1.pbsim/$read2/$read2.RES.txt\e[0m)\n\n");
		my @dedupe_output = `cat $input1.pbsim/$read2/$read2.RES.txt`;
		for (my $j = 0; $j < @dedupe_output; $j++) {
			LOG($outLog, "\t\e[0;32m$dedupe_output[$j]");
		}
		LOG($outLog, "\e[0m\n\n");
	}
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
	LOG($outLog, "\n\e[1;33m====================================== RESULT =======================================\e[0m\n\n\n") if $total_done == 0;
	LOG($outLog, "(also created in .csv format (google sheet friendly) at \e[1;33m$input1.dedupe_stat.csv\e[0m)\n\ndedupe2.sh $dedupe_param2 in=<in> out=<out>\n\n") if $total_done == 0;
	LOG($outLog, "\e[1;33m$total_done\e[0m. line=$linecount, holenum=\e[1;33m$read2\e[0m, seqlen=\e[0;36m$lenseq\e[0m, totalsimreads=\e[1;34m$tot\e[0m, dup=\e[1;35m$dup\%\e[0m, con=\e[1;36m$con\%\e[0m, rem=\e[1;32m$res\%\e[0m\n");
	LOG($outLog, "\e[1;33m$total_done\e[0m. line=$linecount, holenum=\e[1;33m$read2\e[0m, seqlen=\e[0;36m$lenseq\e[0m, maxedits=$maxedits, totalsimreads=\e[1;34m$tot\e[0m, dup=\e[1;35m$dup\%\e[0m, con=\e[1;36m$con\%\e[0m, rem=\e[1;32m$res\%\e[0m\n","NA");
	print $outStat "$total_done,$read2,$lenseq,$maxedits,$tot,$dup,$con,$res\n";
	push(@{$data{tot}}, $tot);
	push(@{$data{dup}}, $dup);
	push(@{$data{con}}, $con);
	push(@{$data{res}}, $res);
	push(@{$data{lenseq}}, $lenseq);
	$total_done ++;
#	DIELOG($outLog, "RES: tot=$tot dup=$dup con=$con res=$res\n") if
#	LOG($outLog, "Done $total_done\n") if $total_done % 10 == 0;
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

LOG($outLog, "\n\e[1;33m================================== STATISTICS OF ALL ===================================\e[0m\n\n\n");
LOG($outLog, "Out of \e[1;33m $total_done \e[0m randomly chosen reads (~ (4 x $simtot) simulated reads each):\n");
LOG($outLog, "length seq  = \e[1;31m$lenseq\e[0m bp +/= \e[1;32m$lenseqsd \e[0m bp\n");
LOG($outLog, "total reads = \e[1;34m$tot\e[0m reads +/= \e[1;32m$totsd\e[0m reads\n");
LOG($outLog, "duplicates  = \e[1;35m$dup\e[0m \% +/= \e[1;32m$dupsd \e[0m \%\n");
LOG($outLog, "containment = \e[1;36m$con\e[0m \% +/= \e[1;32m$consd \e[0m \%\n");
LOG($outLog, "remaining   = \e[1;32m$res\e[0m \% +/= \e[1;32m$ressd \e[0m \%\n");
LOG($outLog, "\n\e[1;33m========================================= DONE ====== ===================================\e[0m\n\n\n");

LOG($outLog, "Output statistics in csv (google sheet friendly): \e[1;33m$input1.dedupe_stat.csv\e[0m\n\n");

draw($res,$ressd);

sub draw {
	my ($num, $numsd) = @_;
	my $orignum = $num;
	my $origsd = $numsd;
my $index = 2;
my $hisd = int(($num+$numsd)/$index+0.5); 
my $losd = int(($num-$numsd)/$index+0.5);
$num = int($num / $index+0.5);
$numsd = int($numsd / $index+0.5);
if ($hisd == $num) {$hisd = $num+1};
if ($losd == $num) {$losd = $num-1};
my $ymax = int(30 / $index+0.5);
my $xmax = int(25 / $index+0.5);
my $xpos1 = int(5 / $index+0.5);
my $xpos2 = int(15 / $index+0.5);
my $xmid = int(($xpos1+$xpos2)/2+0.5);
#LOG($outLog, "$num, $numsd, $ymax, $xmax, $xpos1, $xpos2\n");
LOG($outLog, "\% Remaining\n   \^\n   |\n");
for (my $i = $ymax; $i >= 0; $i--) {
	my $y = $i;
	my $y2 = int(2*$y+0.5);
	#$i == -1 ? LOG($outLog, " ") : 
	($y2 % 5 == 0 and $y2 == 0) ? LOG($outLog, "$y2 -|") : 
	($y2 % 5 == 0 and $y2 < 10) ? LOG($outLog, "$y2 -|") : 
	($y2 % 5 == 0 and $y2 != 0) ? LOG($outLog, "$y2-|") : 
   LOG($outLog, "   |");
#	($i % 5 == 0 and $y == 0) ? LOG($outLog, "$y -|") : 
#	($i % 5 == 0 and $y < 10) ? LOG($outLog, "$y -|") : 
#	($i % 5 == 0 and $y != 0) ? LOG($outLog, "$y-|") : 
#   LOG($outLog, "   |");
	for (my $j = 0; $j <= $xmax; $j++) {
		my $x = $j;
#		(print int($x/10) and next) if $y == -2 and $x % 10 == 0;
#		if ($y == -2 and $x % 10 != 0) {
#			LOG($outLog, "");
#			next;
#		}
		if ($y == int($num/2) and $x == $xpos2+2) {
			LOG($outLog, "$orignum \% +/- $origsd \% Remaining");
		}
		elsif ($x == $xmid and $y >= $losd and $y <= $hisd) {
			LOG($outLog, "|") if $y >= $losd and $y <= $hisd;
		}
		elsif ($x > $xpos1 and $x < $xpos2 and ($y == 0 or $y == $num+1)) {
			LOG($outLog, "_");
		}
		elsif (($x == $xpos1 or $x == $xpos2) and $y >= 0 and $y <= $num) {
			LOG($outLog, "|");
		}
		elsif ($y == 0) {
			LOG($outLog, "_");
			next;
		}
		else {
			LOG($outLog, " ");
		}
	}
	LOG($outLog, "\n");
}
}


__END__
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
   #LOG($outLog, "VALUE @value\nMEAN $mean\n");
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
   #LOG($outLog, "VALUE @value\nMEAN $mean\n");
   for (my $i = 0; $i < @value; $i++) {
      $sd += (($value[$i] - $mean)**2) / (@value-1);
   }
   $sd = sqrt($sd);
   return($sd);
}


__END__
     

Ideally 100\% simulated reads should either be duplicates (100\% identical) or containment (similar enough to be duplicates based on maxedits and minidentity)
However, even with maximum settings, only 85-90\% maximum reads are marked as duplicates/containment.
Therefore the goal is to find parameters that'll call 85-90\% reads as duplicates/containments, 
leaving 10-15\% remaining reads (false negative is at most 10-15\%).

On containment parametes:
- minidentity specifies how many mismatches to still be considered \"similar enough\".
  Default is 98 \%, as expected pacbio ccs mismatches is 2\% at most.
- maxedits specifies how many indels to still be considered \"similar enough\".
  Default -e is 2\% of read length, as expected pacbio ccs indels is 2\% at most.
  2\% is also the threshold which consistently call maximum \% simulated reads into containments
  (leaving 10-15\% \"duplicates\"). 
  More than 2\% will cause no/minuscule change on \% containment in simulated reads, 
  but on the real data it'll significantly call more reads into containment (false positives)
  Less than that it'll sometimes cause significantly less \% containment in simulated reads, hence false negatives.
  1-2\% causes some reads to have only 60-70\% containments from testing.


  \e[1;32mTherefore, we should use maxedits 2\% (more false positives) to be on the more conservative side, making sure
  that our actual data only have at most 10-15\% duplicates (false negatives).\e[0m

-m : [default: 98] percent identity for dedupe.sh minidentity parameter.


	if (defined $opt_p) {
		my @dedupe_param = split(",", $opt_p);
		foreach my $params (@dedupe_param) {
			my ($param, $value) = split("=", $params);
			DIELOG($outLog, "Undefined param=$param or value=$value in params=$params\n") if not defined $param or not defined $value;
			$param = "mo" if $param eq "minoverlap";
			$param = "nam" if $param eq "numaffixmaps";
			$param = "e" if $param eq "maxedits";
			$param = "mid" if $param eq "minidentity";
			next if $value eq "NA";
			if ($param eq "e") {
				$value = $maxedits;
			}
			$param{$param} = "$param=$value";
		}
	}
	my $dedupe_param2 = $dedupe_param;
	$dedupe_param2 
	foreach my $param (sort keys %param) {
		$dedupe_param2 .= " $param=$param{$param}";
	}

