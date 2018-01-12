#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_f $opt_S $opt_c $opt_C $opt_o);
getopts("s:i:g:f:S:cCo:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}

use myFootLib;
my ($folder) = $opt_f;
my ($outDir) = $opt_o;
die "\nusage: $YW$0$N $CY-f [folder of -n footLop.pl]$N $LGN-o$N [output dir]\n\n" unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_f) == 1;
die "\nplease define output (-o)\n" if not defined $opt_o;
makedir($outDir);
my $logFile = "$folder/logFile.txt";

print "Logfile = $LRD$logFile$N\n"; 
my ($samFile, $seqFile, $genez) = parse_logFile($logFile); 
print "Checking sam File =$LCY$samFile$N=\n";
check_file($samFile, "sam"); 
print "Checking seq File =$LCY$seqFile$N=\n";
check_file($seqFile, "seq"); 
my ($folder1, $fileName1) = getFilename($samFile, "folder"); 
my %refs = %{parse_seqFile($seqFile)}; 

foreach my $chr (sort keys %refs) {
	$genez->{$chr} = @{$refs{$chr}};
}
my %out;
my %data; my $cons; my %strand;
my $linecount = 0;
my ($samFolder, $samName) = getFilename($samFile, "folderfull");
my $debugFile = "$outDir/debug.txt";
my $outSam = "$outDir/$samName.fixed";
open (my $outsam, ">", "$outSam") or die "Cannot write to $outSam: $!\n";
open (my $outdebug, ">", "$debugFile") or die "Cannot write to $debugFile: $!\n";
my ($total_read) = `awk '\$2 == 0|| \$2 == 16 {print}' $samFile | wc -l` =~ /^(\d+)$/;
$linecount = 0;
open (my $in1, "<", $samFile) or die "Cannot read from $samFile: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	$linecount ++;
	print date() . "\t$0: Parsed $LGN$linecount$N / $LCY$total_read$N\n" if $linecount % 50 == 0;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted, @others) = @arr;
	$chr = uc($chr);
	my $others = join("\t", @others); $others = @others == 0 ? "" : "\t$others";
	my @ref1 = defined $refs{$chr} ? @{$refs{$chr}} : die "Can't find gene $chr in $seqFile!\n";
	my ($ref2, $seq2, $poz, $seqborder0, $seqborder1) = parse_samFile($line, \@ref1);
	my %poz = %{$poz};
	my @seq1 = split("", $seqs);
	my %bad = %{get_bad_region($ref2, $seq2, $seqborder0, $seqborder1)};
	my ($ref3, $seq3, $bad3);
	for (my $i = 0; $i < @{$ref2}; $i++) {
		my $bad2 = ($i < $seqborder0 or $i >= $seqborder1) ? " " : defined $bad{$i} ? $LRD . "!" . $N : " ";
		if ($ref2->[$i] ne "-") {
			push(@{$ref3}, $ref2->[$i]);
			push(@{$seq3}, $seq2->[$i]);
			push(@{$bad3}, $bad2);
		}
	}
	my	($CTcons, $CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1) = det_C_type($ref3, $seq3, $bad3, $seqborder0, $seqborder1);
	my ($refPrint, $seqPrint) = colorconv($ref3, $seq3);
	my $CTPrint = join("", @{$CTcons});
	my $newstrand = $CT1 > $GA1 ? 0 : $GA1 > $CT1 ? 16 : $strand;
	my $type = "3_NONE" if ($CT1 <= 5 and $GA1 <= 5);
	$type = "6_BOTH" if $CT1 == $GA1 or ($CT1 > 5 and $GA1 > 5 and $GA1 >= $CT1 * 0.9 and $GA1 <= $CT1 * 1.1 and $CT1 >= $GA1 * 0.9 and $CT1 <= $GA1 * 1.1);
	$type = "2_WNEG" if $CT1 < $GA1;
	$type = "4_WPOS" if $CT1 > $GA1;
	$type = "5_SPOS" if ($CT1 > 15 and $GA1 < 5) or ($GA1 >= 5 and $CT1 / $GA1 > 3) or ($GA1 >= 20 and $CT1 / $GA1 >= 2);
	$type = "1_SNEG" if ($GA1 > 15 and $CT1 < 5) or ($CT1 >= 5 and $GA1 / $CT1 > 3) or ($CT1 >= 20 and $GA1 / $CT1 >= 2);
	print $outdebug ">$read,$type,OldStrand=$strand,NewStrand=$newstrand,$chr,CT0=$CT0,CC0=$CC0,GA0=$GA0,GG0=$GG0,CT1=$CT1,CC1=$CC1,GA1=$GA1,GG1=$GG1\n";
	print $outdebug "$refPrint\n";
	print $outdebug "$seqPrint\n";
	print $outdebug "$CTPrint\n";
	print $outsam "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\n";
#	print "$read\t$chr\tstrand=$strand, new=$newstrand\n" . join("", @{$ref3}) . "
#	CC = $CC1 / $CC0
#	CT = $CT1 / $CT0
#	GG = $GG1 / $GG0
#	GA = $GA1 / $GA0
#	REF: $refPrint
#	SEQ: $seqPrint
#	CON: " . join("", @{$CTcons}) . "\n";
}
close $outsam;

foreach my $strand (sort keys %strand) {
	my @types = ("same","diff");
	print $outdebug "$strand: ";
	foreach my $type (@types[0..1]) {
		my $total = $strand{$strand}{$type}; $total = 0 if not defined $total;
		print $outdebug "$type=$total,";
		my $CT = $strand{$strand}{CT}{$type};
		my ($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $CT) {
			$tmm = int(1000*tmm(@{$CT})+0.5)/1000;
			$mean = int(1000*mean(@{$CT})+0.5)/1000;
			$tmmse = int(1000*tmmse(@{$CT})+0.5)/1000;
			$meanse = int(1000*se(@{$CT})+0.5)/1000;
		}
		print $outdebug "CT=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse, ";
		my $tot = $strand{$strand}{tot}{$type}; 
		($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $tot) {
			$tmm = int(1000*tmm(@{$tot})+0.5)/1000;
			$mean = int(1000*mean(@{$tot})+0.5)/1000;
			$tmmse = int(1000*tmmse(@{$tot})+0.5)/1000;
			$meanse = int(1000*se(@{$tot})+0.5)/1000;
		}
		print $outdebug "tot=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse\n";
	}
}
exit 0;

sub check_file {
	my ($file, $type) = @_;
	die "$type file $file does not exist!\n" if ex($file) == 0;
	die "$type file $file is empty!\n" if -s $file == 0;
	my $filetype = `file -b --mime-type $file`; chomp($filetype);
	my @line = ($file =~ /\.(rmdup|bam)$/ or $filetype =~ /(gzip|binary)/) ? `samtools view $file | head -n 200` : `head -n 200 $file`;
	my $check = 0;
	for (my $i = 0; $i < @line; $i++) {
		my $line = $line[$i];
		chomp($line); my @arr = split("\t", $line);
		if ($type eq "sam") {
			$check = 2 if @arr > 7; last if $check == 2;
		}
		if ($type eq "seq") {
			$check = 1 if $line =~ /^>/;
			last if $i == @line-1;
			$i ++; $line = $line[$i];
			$check = 2 if $line !~ /^>/ and $line =~ /^[ACTGUN]+$/ and $check == 1;
			last if $check == 2;
		}
	}
	die "$file does not look like a $type file!\n" if $check != 2;
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
		#print "REF $def length=" . length($seq) . " SEQ = @seq\n\n" if $def eq "CALM3";
   }
   close $in;
	return(\%ref);
}

sub parse_logFile {
	my ($logFile) = @_;
	die "\n\nCan't find $logFile! Please run footLoop.pl first before running this!\n\n" if not -e $logFile;
	print "LOGFILE=$logFile\n";
	my ($samFile, $seqFile, $genez);
	if (-e $logFile) {
		my @line = `cat $logFile`;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i]; chomp($line);
			if ($line =~ /^!samFile=/) {
				($samFile) = $line =~ /^!samFile=(.+)$/;
			}
			if ($line =~ /^!seqFile=/) {
				($seqFile) = $line =~ /^!seqFile=(.+)$/;
			}
			if ($line =~ /gene=.+length=/) {
				my ($gene, $length) = $line =~ /^.+gene=(.+) length=(\d+)$/;
				die if not defined $gene or not defined $length;
				$genez->{$gene} = $length;
				print "gene $gene = $length bp\n";
			}
#			if ($line =~ /2. Parsing in sequence for genes from sequence file/) {
#				for (my $j = $i+1; $j < @line; $j++) {
#					#print "\t$line\n";
#					last if $line =~ /SUCCESS.+Sequence has been parsed from fasta file/;
#					}
#				}
#			}
			last if $line =~ /footLoop_2_sam_to_peak.pl/;
		}
	}
	return($samFile, $seqFile, $genez);
}

sub get_bad_region {
	my ($ref2, $seq2, $seqborder0, $seqborder1) = @_;
	my %bad;
	my @bad; my $prints;
	for (my $i = $seqborder0; $i < $seqborder1; $i++) {
		my $badcount = 0;
		my $reftemp = "";
		my $seqtemp = "";
		for (my $j = $i; $j < $seqborder1; $j++) {
			die "Undefined ref2 at j=$j\n" if not defined $ref2->[$j];
			$reftemp .= $ref2->[$j];
			die "Undefined seq2 at j=$j\n" if not defined $seq2->[$j];
			$seqtemp .= $seq2->[$j];
			if ($ref2->[$j] eq "-" or $seq2->[$j] eq "-") {
				$badcount ++;
			}
			else {
				last;
			}
		}
		my $beg = $i - $badcount < 0 ? 0 : $i - $badcount;
		my $end = $i + 2*$badcount > @{$ref2} ? @{$ref2} : $i + 2*$badcount;
		$prints .= join("", (" ")x$beg) . join("", ("!")x($end-$beg+1)) . "($beg-$end)\n" if $badcount > 0;
		$prints .= join("", (" ")x$i) . "$reftemp\n" if $badcount > 0;
		$prints .= join("", (" ")x$i) . "$seqtemp\n" if $badcount > 0;
		for (my $j = $beg; $j <= $end; $j++) {
			last if $badcount == 0;
			$bad[$j] = 0;
			$bad{$j} = 1;
		}
		$i += $badcount;
	}
	for (my $i = 0; $i < @bad; $i++) {
		$bad[$i] = " " if not defined $bad[$i];
	}
	return(\%bad);
}

sub parse_samFile {
	my ($line, $refs) = @_;
	my @refs = @{$refs};
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted) = split("\t", $line);
	my @seq = split("",$seqs);
	my ($num, $alp, $lengthseq) = parse_cigar($cigar); die if not defined $num;
	my @num  = @{$num}; my @alp = @{$alp};
	my @ref0 = @refs[0..$pos-2];
	my @ref = @refs[$pos-1..@refs-1];
	my ($seq, $ref, $seqpos, $refpos) = (\@seq, \@ref, 0, 0);
	my %pos;
	my $lengthref = @ref; my $insref = 0; my $insseq = 0;
	for (my $i = 0; $i < @num; $i++) {
		#if ($alp[$i] =~ /^(I|D)$/) {
		#for (my $j = 0; $j < $num[$i]; $j++) {
		($ref) = $alp[$i] eq "I" ? ins($ref, $refpos, "-", $num[$i], "ref") : $ref;
		($seq) = $alp[$i] eq "D" ? ins($seq, $seqpos, "-", $num[$i], "seq") : $seq; 
		$refpos += $num[$i];
		$seqpos += $num[$i];
#				$refpos ++;
#				$seqpos ++;
		$insref += $num[$i] if $alp[$i] eq "I";
		$insseq += $num[$i] if $alp[$i] eq "D";
				#$pos{seq}{$seqpos-$insseq+$pos-1} = $refpos-1+$pos-1;
#				#$pos{ref}{$refpos-$insref+$pos-1} = $seqpos-1+$pos-1;
#			}
#		}
#		else {
#		}
	}
	my $refend = $refpos - $insref;
	@ref = (@ref0, @{$ref});
	my $seqborder0 = $pos-1;
	my $seqborder1 = $pos - 1 + @{$seq};
	@seq = ((("-") x ($pos-1)), @{$seq}, (("-")x($lengthref-$refend)));
	return(\@ref, \@seq, \%pos, $seqborder0, $seqborder1);
}

sub det_C_type {
   my ($ref, $seq, $bad, $seqborder0, $seqborder1) = @_;
   my @ref = @{$ref};
   my @seq = @{$seq};
	my @bad = @{$bad};
   my (@top);
   @ref = ("-","-",@ref,"-","-");
   @seq = ("-","-",@seq,"-","-");
   @bad = (" "," ",@bad," "," ");
   $seqborder0 += 2;
   $seqborder1 += 2;
   my $len = @ref;
	my ($CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1) = (0,0,0,0,0,0,0,0);
   for (my $i = 2; $i < $len-2; $i++) {
      my ($beg, $end, $pos) = ($i-2, $i+2, $i-2); #pos is for top and bot position
      my ($top1, $top2, $top3, $top4, $top5) = @ref[$beg..$end]; #NNCNN
      my ($bot1, $bot2, $bot3, $bot4, $bot5) = @seq[$beg..$end]; #NNGNN
		my $chunkbot = join("", @seq[$beg..$end]);
      if ($i < $seqborder0 or $i >= $seqborder1) {
			$top[$pos] = "-";
      }
      elsif ($top3 eq "C") {
			my $top;
			if ($bad[$i] eq " ") {
	         $top = $top4 eq "G" ? "E" : $top5 eq "G" ? "D" : "C";
         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
			else {
	         $top = $top4 eq "G" ? "3" : $top5 eq "G" ? "2" : "1";
         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
         $top = $bot3 eq "T" ? lc($top) : ($top =~ /^(C|D)$/ and $bot3 eq "-") ? "B" : ($top =~ /^E$/ and $bot3 eq "-") ? "F" : $top;
         $top =~ tr/CDE123/MNOPQR/ if $bot3 !~ /^(C|T)$/; #not CC or CT.
			#CH=B(-)CDcdMN, CG=F(-)EeO
         $top[$pos] = $top;
      }
      elsif ($top3 eq "G") {
			my $top;
			if ($bad[$i] eq " ") {
	         $top = $top2 eq "C" ? "I" : $top1 eq "C" ? "H" : "G";
	         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
			else {
	         $top = $top2 eq "C" ? "6" : $top1 eq "C" ? "5" : "4";
	         #             --- ..CG. ? --- : ---- ..CHG ? --- : ..CHH
			}
         $top = $bot3 eq "A" ? lc($top) : ($top =~ /^(G|H)$/ and $bot3 eq "-") ? "J" : ($top eq "I" and $bot3 eq "-") ? "K" : $top;
         $top =~ tr/GHI456/UVWXYZ/ if $bot3 !~ /^(G|A)$/; # not GG or GA. 
			#GH=J(-)GHghUV, GC=K(-)IiW
			$top[$pos] = $top;
      }
      else {
			$top[$pos] = $bot3 eq "-" ? "_" : "$bot3";
      }
		$CC0 ++ if "$top3$bot3" eq "CC"; 
		$GG0 ++ if "$top3$bot3" eq "GG"; 
		$CT0 ++ if "$top3$bot3" eq "CT"; 
		$GA0 ++ if "$top3$bot3" eq "GA";
		if ($bad[$i] eq " ") {
			my $left = $i - 5 < 0 ? 0 : $i - 5;
			my $rite = $i + 5 >= @ref ? @ref - 1: $i + 5;
			my $chunk = join("", @seq[$left..$rite]);
			if ($chunk !~ /\-/) {
				$CT1 ++ if "$top3$bot3" eq "CT";
				$GA1 ++ if "$top3$bot3" eq "GA";
				$CC1 ++ if "$top3$bot3" eq "CC"; 
				$GG1 ++ if "$top3$bot3" eq "GG"; 
			}
		}
   }
   return(\@top, $CC0, $GG0, $CC1, $GG1, $CT0, $GA0, $CT1, $GA1);
}
sub ins {
	my ($arr, $pos, $ins, $total, $type) = @_;
	my @ins = ($ins) x $total;
	my @arr = @{$arr};
	my @arr0 = @arr[0..$pos-1];
	my @arr1 = @arr[$pos..@arr-1];
	die if @arr0 == 0;
	@arr = (@arr0, @ins, @arr1);
	return(\@arr);
}

sub getConv {
	my ($converted) = @_;
	($converted) = $converted =~ /^.+:([\.A-Z]+)$/i; my $length2 = length($converted);
	my ($X) = $converted =~ tr/X/X/;
	my ($U) = $converted =~ tr/U/U/;
	my ($H) = $converted =~ tr/H/H/;
	my ($Z) = $converted =~ tr/Z/Z/;
	my ($x) = $converted =~ tr/x/x/;
	my ($u) = $converted =~ tr/u/u/;
	my ($h) = $converted =~ tr/h/h/;
	my ($z) = $converted =~ tr/z/z/;
	my ($dot) = $converted =~ tr/\./\./;
	my $length = ($X+$H+$Z+$U+$x+$h+$z+$u+$dot);
	die "X=$X, H=$H, Z=$Z, U=$u, x=$x, h=$h, z=$z, u=$u, dot=$dot, length=$length, length2=$length2\n\n" if $length ne $length2;
	$data{stat} = "X=$X, H=$H, Z=$Z, U=$u, x=$x, h=$h, z=$z, u=$u, dot=$dot, lengthsum=$length, lengthcol14=$length2";
	my $conv = $x + $h + $z + $u;
	my $notc = $X + $H + $Z + $U;
	my $nonc = $dot;
	my @converted = split("", $converted);
	return(\@converted);

}
sub check_chr_in_sam {
	my ($samFile) = @_;
	open (my $in2, "cut -f2,3 $samFile|") or die "Cannot read from $samFile: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		my @arr = split("\t", $line);
		next if @arr == 0;
		next if $arr[0] !~ /^\d+$/;
		$linecount ++;
		my ($strand, $chr) = @arr; $chr = uc($chr);
		next if $strand eq 4;
		die "Can't find gene $chr in $seqFile!\n" if not defined $refs{$chr};
		next;
	}
}
__END__
sub colorseq {
   my ($seq) = @_;
   my @seq = $seq =~ /ARRAY/ ? @{$seq} : split("", $seq);
	my $seq3 = "";
   foreach my $seq2 (@seq[0..@seq-1]) {
		$seq3 .= color($seq2);
   }
	return ($seq3);
}
sub colorseqCG {
   my ($seq) = @_;
   my @seq = $seq =~ /ARRAY/ ? @{$seq} : split("", $seq);
	my $seq3 = "";
   foreach my $seq2 (@seq[0..@seq-1]) {
		$seq3 .= colorCG($seq2);
   }
	return ($seq3);
}

sub color {
   my ($nuc) = @_;
   $nuc = $nuc eq "A" ? "$LRD$nuc$N" : $nuc eq "G" ? "$LGN$nuc$N" : $nuc eq "C" ? "$YW$nuc$N" : $nuc eq "T" ? "$LCY$nuc$N" : "$LPR$nuc$N";
   return $nuc;
}

sub colorCG {
   my ($nuc) = @_;
   $nuc = $nuc eq "C" ? "$LGN$nuc$N" : $nuc eq "G" ? "$LGN$nuc$N" : "$GR$nuc$N";
   return $nuc;
}

sub colorconv {
   my ($seq1, $seq2) = @_;
	my @seq1 = $seq1 =~ /ARRAY/ ? @{$seq1} : split("", $seq1);
	my @seq2 = $seq2 =~ /ARRAY/ ? @{$seq2} : split("", $seq2);
	my ($res1, $res2);
	for (my $i = 0; $i < @seq1; $i++) {
		my $dinuc = $seq1[$i] . $seq2[$i];
   	$res1 .= $dinuc eq "CT" ? "${LRD}C$N" : $dinuc eq "GA" ? "${LRD}G$N" : ($seq1[$i] =~ /^(C|G)$/ and $seq1[$i] ne $seq2[$i]) ? "${YW}$seq1[$i]$N" : colorCG($seq1[$i]);
   	$res2 .= $dinuc eq "CT" ? "${LPR}T$N" : $dinuc eq "GA" ? "${LPR}A$N" : ($seq2[$i] =~ /^(C|G)$/ and $seq2[$i] ne $seq1[$i]) ? "${YW}$seq2[$i]$N" : colorCG($seq2[$i]);
	}
	return($res1, $res2);
}

__END__

