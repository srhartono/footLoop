#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_t $opt_l $opt_i $opt_c $opt_m);
getopts("vt:l:i:cm:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";

my ($window, $thres, $input1, $minCG) = ($opt_l, $opt_t, $opt_i, $opt_m);
die "\nUsage: $YW$0$N $LGN-l$N <length>$LGN -t$N <\% threshold> $CY-i $N<folder from -n footLoop.pl>

-i <output folder from -n parameter in footLoop.pl>
-t [50]  <\% threshold>
-l [100] <window size bp>
-c <toggle include CpG in calculation>
-m <min C or G for peak>

" unless defined $opt_i and -d $opt_i;

my %opts = ("i" => $opt_i, "t" => $opt_t, "l" => $opt_l, "c" => $opt_c, "m" => $opt_m);
$window = 100 if not defined $window;
$thres  = 50  if not defined $thres;
$minCG  = 5   if not defined $minCG;

my $logFile = "$opt_i/logFile.txt";
print date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$logFile$N from footLoop.pl does not exist!: $!\n" and exit 1 if not -e $logFile;
my ($samFile, $seqFile, $genez) = parse_logFile($logFile);
print date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$seqFile$N from footLoop.pl does not exist!: $!\n" and exit 1 if not defined $seqFile or not -e $seqFile;
my $refs = parse_seqFile($seqFile);
print date() . ":$LRD FATAL ERROR$N: $YW$0$N Cannot get reference sequence from $seqFile!: $!\n" and exit 1 if not defined $refs;
if (-d $opt_i) {
	open (my $outLog, ">", "$opt_i/logFile2.txt") or print STDERR date() . ":$LRD FATAL ERROR$N: $YW$0$N Failed to write into $CY$opt_i/logFile2.txt$N: $!\n" and exit 1;
	record_options(\%opts, $outLog);
	my @inputs = <$opt_i/0_orig/*.orig>;
	print $outLog date() . ":$LRD FATAL ERROR$N: $YW$0$N There is no orig files in $CY$opt_i/0_orig/$N: $!\n" and exit 1 if @inputs == 0;
	for (my $i = 0; $i < @inputs; $i++) {
		next if $inputs[$i] !~ /CALM3_/;
		print $outLog date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$inputs[$i]$N is empty or does not exist!: $!\n" and next if -s $inputs[$i] == 0;
		my $status = parse($inputs[$i], $outLog, $i, $refs); $status = "FAIL" if not defined $status;
		print "$inputs[$i] Status = $status\n";
		print "\n\n";# and exit 1;
	}
}
elsif (-e $opt_i) {
	die "-i must be folder output of footLoop.pl -n!\n";
	open (my $outLog, ">", "$opt_i.logFile2.txt") or print STDERR date() . ":$LRD FATAL ERROR$N: $YW$0$N Failed to write into $CY$opt_i.logFile2.txt$N: $!\n" and exit 1;
	record_options(\%opts, $outLog);
	my $status = parse($opt_i, $outLog, 1);
	print "Status = $status\n";
}

sub parse {
	my ($input1, $outLog, $fileNum, $refs) = @_;
	($input1) = getFullpath($input1);
	my ($folder1, $fileName1) = getFilename($input1, "folder");
	my ($gene, $strand) = $fileName1 =~ /^(.+)_(Pos|Neg).orig$/;
	print $outLog date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$input1$N gene or strand can't be extracted from $fileName1!: $!\n" and return 1 if not defined $gene or not defined $strand;
	my $ref = $refs->{$gene};
	print $outLog date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$input1$N gene=$gene cannot get sequence from $seqFile!: $!\n" and return 2 if not defined $ref;

	print "\t$fileNum. Processing $LGN$fileName1$N ($LCY$input1$N)\n\n";
	my %peak;
	my $lineCount = 0;
#	print "Here ";# and return "${LGN}SUCCESS$N";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	open (my $out1, ">", "$input1.Rdata") or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		$lineCount ++;
		my ($read, $type, $seq) = split("\t", $line);
		print $outLog date() . ":$LRD FATAL ERROR$N: $YW$0$N $LCY$input1$N lineCount=$lineCount undefined! read=$read type=$type seq=$seq: $!\n" and return 3 if not defined $read or not defined $type or not defined $seq;
		
		my @seq = split("", $seq);
		last if $lineCount > 2000;
		print "$input1: Done $lineCount\n" if $lineCount % 250 == 0;
		my %count; my $printR;
		for (my $i = 0; $i < @seq-$window; $i+= int($window/4+0.5)) {
			my $seqChunk = join("", @seq[$i..($i+$window-1)]);
			my $refChunk = join("", @{$ref}[$i..($i+$window-1)]);
#			($count{C}{tot}) = $seqChunk =~ tr/Cc/Cc/;
#			($count{G}{tot}) = $seqChunk =~ tr/Gg/Gg/;
#			($count{C}{use}) = $seqChunk =~ tr/c/c/; # all C except those at deleted region
#			($count{G}{use}) = $seqChunk =~ tr/g/g/; # all C except those at deleted region
			($count{CH}{tot})    = $seqChunk =~ tr/CDcdMNB/CDcdMNB/; # all CH
			($count{CH}{use})    = $seqChunk =~ tr/cd/cd/; # conv CH
			($count{CG}{tot})    = $seqChunk =~ tr/EeOF/EeOF/; # all CG
			($count{CG}{use})    = $seqChunk =~ tr/e/e/; # conv CG
			($count{GH}{tot})    = $seqChunk =~ tr/GHghUVJ/GHghUVJ/; # all HG
			($count{GH}{use})    = $seqChunk =~ tr/gh/gh/; # conv HG
			($count{GC}{tot})    = $seqChunk =~ tr/IiWK/IiWK/; # all HG
			($count{GC}{use})    = $seqChunk =~ tr/i/i/; # conv HG
#			($count{CH}{tot})    = $seqChunk =~ tr/Cc/Cc/; # all CH
#			($count{CH}{use})    = $seqChunk =~ tr/c/c/; # conv CH
#			($count{GH}{tot})    = $seqChunk =~ tr/Gg/Gg/; # all HG
#			($count{GH}{use})    = $seqChunk =~ tr/g/g/; # conv HG
=comment
			while ($seqChunk =~ /(MZ|ZZ|ZN|[A-LO-Y]Z$|^Z[A-LO-Y])/ig) {#([A-Z\-\. ]Z[A-Z\-\. ])/ig) {
				my $curr = $&; my @curr = split("", $curr);
				my $C = $curr =~ /z./ ? "c" : $curr =~ /Z./ ? "C" : $curr[0];
				my $G = $curr =~ /.z/ ? "g" : $curr =~ /.Z/ ? "G" : $curr[1];
				my $type = "$C,$G";
				($count{CG}{tot}) ++ if $C =~ /[CcMm]/;
				($count{CH}{tot}) -- if $C =~ /[Mm]/;
				($count{GC}{tot}) ++ if $G =~ /[GgNn]/;
				($count{GH}{tot}) -- if $G =~ /[Nn]/;
				($count{CG}{use}) ++ if $C =~ /[cm]/;
				($count{GC}{use}) ++ if $G =~ /[gn]/;
				($count{CH}{use})  -- if $C =~ /[m]/;
				($count{GH}{use})  -- if $G =~ /[n]/;
				#print "$curr: $type\n";
			}
=cut
			my @keys = qw(tot use);
			my @nucs = qw(CG GC CH GH);#C G);
			foreach my $nuc (@nucs[0..@nucs-1]) {
				foreach my $key (@keys[0..@keys-1]) {
					$count{$nuc}{$key} = 0 if not defined $count{$nuc}{$key};
				}
			}
#			undef $opt_c;
			my $Ctot = $count{CH}{tot};
			my $Gtot = $count{GH}{tot};
#			my $Ctot = defined $opt_c ? $count{CG}{tot} + $count{CH}{tot} : $count{CH}{tot};
#			my $Gtot = defined $opt_c ? $count{GC}{tot} + $count{GH}{tot} : $count{GH}{tot};
			my $conC = $Ctot < $minCG ? 0 : int(0.5+1000*($count{CH}{use}) / ($count{CH}{tot}))/10;
			my $conG = $Gtot < $minCG ? 0 : int(0.5+1000*($count{GH}{use}) / ($count{GH}{tot}))/10;
#			my $conC = $Ctot < $minCG ? 0 : defined $opt_c ? int(0.5+1000*($count{CG}{use} + $count{CH}{use}) / ($count{CG}{tot} + $count{CH}{tot}))/10 : int(0.5+1000*($count{CH}{use}) / ($count{CH}{tot}))/10;
#			my $conG = $Gtot < $minCG ? 0 : defined $opt_c ? int(0.5+1000*($count{GC}{use} + $count{GH}{use}) / ($count{GC}{tot} + $count{GH}{tot}))/10 : int(0.5+1000*($count{GH}{use}) / ($count{GH}{tot}))/10;
			$printR->{CHper} .= "\t$conC";
			$printR->{CHtot} .= "\t$Ctot";
			$printR->{GHper} .= "\t$conG";
			$printR->{GHtot} .= "\t$Gtot";
#			print $out1 "$read\tCH=$conC \% of $Ctot\tGH=$conG % of $Gtot";
			$peak{$read}{C}{$i} = $conC;
			$peak{$read}{G}{$i} = $conG;
			$Ctot = $count{CG}{tot} + $count{CH}{tot};
			$Gtot = $count{GC}{tot} + $count{GH}{tot};
#			die "tot C=$Ctot, CG=$count{CG}{tot}, CH=$count{CH}{tot}\n" if $count{CG}{tot} + $count{CH}{tot} == 0 and $Ctot > 0;
			$conC = $Ctot < $minCG ? 0 : int(0.5+1000*($count{CG}{use} + $count{CH}{use}) / ($count{CG}{tot} + $count{CH}{tot}))/10;
#			die "tot G=$Gtot, GC=$count{GC}{tot}, GH=$count{GH}{tot}\n" if $count{GC}{tot} + $count{GH}{tot} == 0 and $Gtot > 0;
			$conG = $Gtot < $minCG ? 0 : int(0.5+1000*($count{GC}{use} + $count{GH}{use}) / ($count{GC}{tot} + $count{GH}{tot}))/10;
			$printR->{CGper} .= "\t$conC";
			$printR->{CGtot} .= "\t$Ctot";
			$printR->{GCper} .= "\t$conG";
			$printR->{GCtot} .= "\t$Gtot";
#			print $out1 "\tCG=$conC \% of $Ctot\tGC=$conG % of $Gtot\n";
#			print "$i $strand $conC $conG\n" if $i == 2238;
#			foreach my $nuc (@nucs[0..@nucs-1]) {
#				my $tot = $count{$nuc}{tot};
#				my $use = $count{$nuc}{use};
#				print "$nuc\t$use/$tot\n";
#			}
#			printChunk($refChunk, $seqChunk) if $i == 2238;
#			last if $i > 2000;# and die;# "died at fileNum=$fileNum file=$input1 i=$i read=$read type=$type seq=\n\n$seq\n\n" if $count{CpG} % 2 != 0;
			
		}
		print $out1 "$read\tCH\tpercent$printR->{CHper}\n" if defined $printR->{CHper};
		print $out1 "$read\tCH\ttotal$printR->{CHtot}\n" if defined $printR->{CHtot};
		print $out1 "$read\tGH\tpercent$printR->{GHper}\n" if defined $printR->{GHper};
		print $out1 "$read\tGH\ttotal$printR->{GHtot}\n" if defined $printR->{GHtot};
		print $out1 "$read\tCG\tpercent$printR->{CGper}\n" if defined $printR->{CGper};
		print $out1 "$read\tCG\ttotal$printR->{CGtot}\n" if defined $printR->{CGtot};
		print $out1 "$read\tGC\tpercent$printR->{GCper}\n" if defined $printR->{GCper};
		print $out1 "$read\tGC\ttotal$printR->{GCtot}\n" if defined $printR->{GCtot};
	}
	my %count;
	close $in1; $count{G}{peak} = 0; $count{C}{peak} = 0;
	my %print; $print{C} = ""; $print{G} .= "";
	foreach my $read (sort keys %peak) {
		$count{tot} ++;
		my $check = 0; my $prints = "";
		if (defined $peak{$read}{C}) {
			foreach my $pos (sort {$a <=> $b} keys %{$peak{$read}{C}}) {
				$prints .= ";$pos,$peak{$read}{C}{$pos}" if $peak{$read}{C}{$pos} > $thres;
				$check = 1 if $peak{$read}{C}{$pos} > $thres;
			}
		}
		$print{C} .= "${YW}C$N: $prints\n" if $check == 1;
		$count{C}{peak} += $check;
		$check = 0; $prints = "";
		if (defined $peak{$read}{G}) {
			foreach my $pos (sort {$a <=> $b} keys %{$peak{$read}{G}}) {
				$prints .= ";$pos,$peak{$read}{G}{$pos}" if $peak{$read}{G}{$pos} > $thres;
				$check = 1 if $peak{$read}{G}{$pos} > $thres;
			}
		}
		$print{G} .= "${LGN}G$N: $prints\n" if $check == 1;
		$count{G}{peak} += $check;
	}
#	print "$print{C}\n" if $print{C} =~ /;/;
#	print "$print{G}\n" if $print{G} =~ /;/;
	return "${LGN}SUCCESS$N: total=$count{tot}, C=$count{C}{peak}, G=$count{G}{peak}";
}


#open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
#close $out1;


sub record_options {
   my ($opts, $outLog) = @_;
   my $optPrint = "$0";
   foreach my $opt (sort keys %{$opts}) {
      next if not defined $opts->{$opt};
      my $val = defined $opts->{$opt} ? $opts->{$opt} : "";
		$val = -e $val ? getFullpath($val) : $val;
      $optPrint .= " -$opt $val";
   }
	my $date = getDate();
	my $uuid = getuuid();
	my $param = "
${YW}Initializing...$N
   
Date       : $date
Run ID     : $uuid
Run script : $optPrint

";

	print $outLog $param;
	print STDERR $param;
}

sub parse_logFile {
   my ($logFile) = @_;
   die "\n\nPlease run footLoop.pl first before running this!\n\n" if not -e $logFile;
   my ($samFile, $seqFile, $genez);
   if (-e $logFile) {
      my @line = `cat $logFile`;
      for (my $i = 0; $i < @line; $i++) {
         my $line = $line[$i]; chomp($line);
         if ($line =~ /^!sam=/) {
            ($samFile) = $line =~ /^!sam=(.+)$/;
         }
         if ($line =~ /^!seq=/) {
            ($seqFile) = $line =~ /^!seq=(.+)$/;
         }
#        print "$line\n" if $i > 20 and $i < 50;
         if ($line =~ /2. Parsing in sequence for genes from sequence file/) {
            for (my $j = $i+1; $j < @line; $j++) {
               #print "\t$line\n";
               last if $line[$j] =~ /SUCCESS.+Sequence has been parsed from fasta file/;
               my ($gene, $length) = $line[$j] =~ /^.+genez=(.+) \(.+\) Length=(\d+)$/;
               die if not defined $gene or not defined $length;
               $genez->{$gene} = $length;
   #           print "parse_logFile: parsed gene=$LCY$gene$N=\n";
            }
         }
         #($threshold) = $line =~ /^\d+\.[ ]+\-t .+Threshold.+[ ]:(\-?\d+\.?\d*)$/ if $line =~ /^\d+\.[ $
         last if ex([$samFile, $seqFile]) == 1 and defined $genez;
         last if $line =~ /SUCCESS.+Sequence has been parsed from fasta file/;
      }
   }
#  if (not -e $logFile or not defined $samFile) {
#     my @sam = <$folder/*.[sb]am>;
#     die "Cannot find a .sam..bam file in folder $CY$folder$N\n\n" if @sam == 0;
#     die "More than one .sam/.bam file in folder $CY$folder$N\n-" . join("\n-", @sam) . "\n\n" if @sam $
#     $samFile = $sam[0] if not -e $logFile;
#     $samFile = $sam[0] if not defined $samFile;
#  }
   #die "Cant find threshold from $logFile!\n" if not defined $threshold;
   return($samFile, $seqFile, $genez);
}

sub parse_seqFile {
   my ($seqFile) = @_;
   open (my $in, "<", $seqFile) or die "Failed to open seqFile $CY$seqFile$N: $!\n";
   my %ref;
   my $fasta = new FAlite($in);
   while (my $entry = $fasta->nextEntry()) {
      my $def = $entry->def; $def =~ s/>//; $def = uc($def);
      my $seq = $entry->seq;
      my @seq = split("", $seq);
      @{$ref{$def}} = @seq;
      #print "REF $def SEQ = @seq\n\n" if $def eq "CALM3";
   }
   close $in;
   return(\%ref);
}

sub printChunk {
	my @seqs = @_;
	my @print;
#	foreach my $seqs (@seqs[0..@seqs-1]) {
#		my @temp = split("", $seqs);
	my $seqCount = 0;
	foreach my $seq (@seqs[0..@seqs-1]) {
		$seqCount ++;
		my $count = 0;
		my @seq = split("", $seq);
		@seq = @{$seq} if $seq =~ /^ARRAY/;
		my $max = @seq > 200 ? 200 : @seq;
		for (my $i = 0; $i < $max; $i += 50) {#@seq-50; $i += 50) {
			my $beg = $i;
			my $end = $i+50 >= $max ? $max : $i + 50;#@seq ? @seq : $i + 50;
			my $seq2 = join("", @seq[$beg..$end-1]);
			$print[$count] .= "$YW $count $N" .  colorseq($seq2) . "\n";
			$count ++;
		}
	}
	for (my $i = 0; $i < @print; $i++) {
		print "$print[$i]\n";
	}
}

sub colorseq {
   my ($seq) = @_;
   my @seq = split("", $seq);
	my $seq3 = "";
   foreach my $seq2 (@seq[0..@seq-1]) {
		$seq3 .= color($seq2);
   }
	return $seq3;
}

sub color {
   my ($nuc) = @_;
   $nuc = $nuc eq "A" ? "$LRD$nuc$N" : $nuc eq "G" ? "$LGN$nuc$N" : $nuc eq "C" ? "$YW$nuc$N" : $nuc eq "T" ? "$LCY$nuc$N" : "$LPR$nuc$N";
   return $nuc;
}


__END__
#		if ($curr =~ /[Zz][ZzNn]/) {
#			$count{CG}{use} ++ if $curr =~ /z[ZzNn]/;
##			$type .= $curr =~ /z[ZzNn]/ ? "c" : "C";
#		}
#		if ($curr =~ /[ZzMm][Zz]/) {
#			$count{GC}{use} ++ if $curr =~ /[ZzMm]z/;
#			$type .= $curr =~ /[ZzMm]z/ ? ",g" : ",G";
#		}
#		my ($C, $G) = split(",", $type);
#		$C = "M" if not defined

