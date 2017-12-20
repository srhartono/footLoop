#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_f $opt_S $opt_c $opt_C $opt_t);
getopts("s:i:g:f:S:cCt:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}

use myFootLib;
my ($folder) = $opt_f;
die "\nusage: $YW$0$N $LGN-t threshold [0-100]$N $CY-f [folder of -n footLop.pl]$N\n\n" unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_f) == 1;
my $logFile = "$folder/logFile.txt";

print "Logfile = $LRD$logFile$N\n"; my ($samFile, $seqFile, $genez, $threshold) = parse_logFile($logFile); $threshold = 100 * $threshold if $threshold < 1; $samFile = $opt_s if (ex($opt_S) == 1); $seqFile = $opt_s if (ex($opt_s) == 1); 
$threshold = $opt_t if defined $opt_t;
check_file($samFile, "sam"); check_file($seqFile, "seq"); my ($folder1, $fileName1) = mitochy::getFilename($samFile, "folder"); 
my %refs = %{parse_seqFile($seqFile)}; 
print "Here\n";
foreach my $chr (sort keys %refs) {
	#print "GENES_ $chr\n";
	$genez->{$chr} = @{$refs{$chr}};
}
print "here2\n";
#155 3M 27I 3M ... 4I 1M 6I 4M:   16 204 GAT-----------GCGAGTAGAGCAGTCGAACATGAGCTGACTCAGGTCACCGA     GCTACGATGTGATGCTTGCACAAGTGATCCA 
#324 6M 10I 3M ... 3M 2I 1M 2I 4M: 0 204  GCAGTC--------GAACATGTAGCTGACTCAGGTCACCGATGTACGGGCCAGAT    GCTACGATGTGATGCTTGCACAAGTGATCCA
#XR:Z:CT XG:Z:GA

my %out;
foreach my $gene (sort keys %{$genez}) {
	my $outPos = "$folder/$gene\_Pos$threshold.txt";
	my $outNeg = "$folder/$gene\_Neg$threshold.txt";
	my $outPosCG = "$folder/$gene\_Pos$threshold\_CG.txt";
	my $outNegCG = "$folder/$gene\_Neg$threshold\_CG.txt";
	$out{$gene}{OUTPOS} = $outPos;
	$out{$gene}{OUTNEG} = $outNeg;
	$out{$gene}{OUTPOSCG} = $outPosCG;
	$out{$gene}{OUTNEGCG} = $outNegCG;
	#open (my $outz, ">", $out{$gene}{OUTPOS}) or die "Cannot write to $out{$gene}{OUTPOS}: $!\n";
	#open (my $outz2, ">", $out{$gene}{OUTNEG}) or die "Cannot write to $out{$gene}{OUTNEG}: $!\n";
	#open (my $outz3, ">", $out{$gene}{OUTPOSCG}) or die "Cannot write to $out{$gene}{OUTPOSCG}: $!\n";
	#open (my $outz4, ">", $out{$gene}{OUTNEGCG}) or die "Cannot write to $out{$gene}{OUTNEGCG}: $!\n";
	#close $outz; close $outz2; close $outz3; close $outz4;
	print "$gene: $outPos\n";
}
my %data; my $cons; my %strand;
my $linecount = 0; my $REF;
my ($samFolder, $samName) = getFilename($samFile, "folderfull");
my $debugFile = "$folder/debug.txt";
my $outSam = "$samName.fixed";
open (my $outsam, ">", $outSam) or die "Cannot write to $outSam: $!\n";
open (my $outdebug, ">", $debugFile) or die "Cannot write to $debugFile: $!\n";
print "Processing $samFile\n";
open (my $in1, "<", $samFile) or die "Cannot read from $samFile: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	$linecount ++;
	print "Parsed $linecount\n" if $linecount % 1000 == 0;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted, @others) = split("\t", $line);
	my $others = join("\t", @others); $others = @others == 0 ? "" : "\t$others";
	my @ref1 = defined $refs{uc($chr)} ? @{$refs{uc($chr)}} : die "Can't find gene $chr in $seqFile!\n";
	my ($ref2, $seq2, $con2, $myc2, $poz, $seqborder0, $seqborder1) = parse($line, \@ref1);
	$chr = uc($chr);
	#print "gene=$chr\n";
	$REF = $ref2 if not defined $REF;
	my %poz = %{$poz};
	my @seq1 = split("", $seqs);
	if (defined $opt_C) {
		last if $linecount > 50; 
		($cons) = make_consensus($ref2, $seq2, $seqborder0, $seqborder1, $cons);
#		next;
	}

	my %bad = %{get_bad_region($ref2, $seq2, $seqborder0, $seqborder1)};
	my %count; my $max = @{$ref2};
	my $bad2;
	for (my $i = 0; $i < @{$ref2}; $i++) {
#	for (my $i = $seqborder0; $i < $seqborder1; $i++) {
		$bad2->[$i] = ($i < $seqborder0 or $i >= $seqborder1) ? " " : defined $bad{$i} ? $LRD . "!" . $N : " ";
	}
	my $VAL; my $VALCG; my ($CT, $GA) = (0,0);
	my ($ref3, $seq3, $con3, $myc3, $bad3) = (); my $tempstrand = $strand;
	for (my $i = 0; $i < @{$ref2}; $i++) {
		if ($ref2->[$i] ne "-") {
			my $valcg = ($bad2->[$i] ne " " and $myc2->[$i] ne " ") ? "B" : ($seq2->[$i] eq "-" and $myc2->[$i] eq " ") ? "." : $myc2->[$i];
			my $val = $valcg;
			$valcg = $valcg eq "." ? 6 : $valcg eq "Z" ? 3 : $valcg eq "z" ? 9 : $valcg =~ /[XHU]/ ? 0 : $valcg =~ /[xhu]/ ? 9 : ($valcg ne "B" and $seq2->[$i] ne "-") ? 2 : $valcg;
			$val = $val eq "." ? 6 : $val =~ /[XHZU]/ ? 0 : $val =~ /[xhzu]/ ? 9 : ($seq2->[$i] ne "-" and $val ne "B") ? 2 : $val;
			$CT++ if $ref2->[$i] eq "C" and $seq2->[$i] eq "T" and $bad2->[$i] eq " ";
			$GA++ if $ref2->[$i] eq "G" and $seq2->[$i] eq "A" and $bad2->[$i] eq " ";
			$tempstrand = $CT > $GA ? 0 : $CT < $GA ? 16 : $strand;
			push(@{$ref3}, $ref2->[$i]);
			push(@{$seq3}, $seq2->[$i]);
			push(@{$myc3}, $myc2->[$i]);
			push(@{$con3}, $con2->[$i]);
			push(@{$bad3}, $bad2->[$i]);
			$val = $val !~ /^[0-9]$/ ? -1 : $val;
			$valcg = $valcg !~ /^[0-9]$/ ? -1 : $valcg;
			push(@{$VAL}, $val);
			push(@{$VALCG}, $valcg);
		}
	}
	#print "BEFORE: CT=$CT, GA=$GA\n" if $tempstrand ne $strand;
	my ($VAL2, $VALCG2);
	my	($cons_CT) = det_C_type($ref3, $seq3, $con3, 0, $seqborder0, $seqborder1);
	my	($cons_GA) = det_C_type($ref3, $seq3, $con3, 16, $seqborder0, $seqborder1);
	
	next;
###################### END ######################

	my	($myc4) = det_C_type($ref3, $seq3, $con3, $tempstrand, $seqborder0, $seqborder1);
	if ($tempstrand ne $strand) {
		($CT, $GA) = (0,0);
		my $REFTEMP = $ref3; my $SEQTEMP = $seq3; my $BADTEMP = $bad3;
		my ($ref3, $seq3, $con3, $myc3, $bad3, $VAL, $VALCG) = (); undef $VAL; undef $VALCG;
		die if defined $VAL or (defined $VAL and defined $VAL->[0]);
		for (my $i = 0; $i < @{$REFTEMP}; $i++) {
			if ($REFTEMP->[$i] ne "-") {
				my $valcg = $myc4->[$i];
				$valcg = ($BADTEMP->[$i] ne " " and $myc4->[$i] ne " ") ? "B" : ($SEQTEMP->[$i] eq "-" and $myc4->[$i] eq " ") ? "." : $myc4->[$i];
				my $val = $valcg;
				$valcg = $valcg eq "." ? 6 : $valcg eq "Z" ? 3 : $valcg eq "z" ? 9 : $valcg =~ /[XHU]/ ? 0 : $valcg =~ /[xhu]/ ? 9 : ($valcg ne "B" and $SEQTEMP->[$i] ne "-") ? 2 : $valcg;
				$val = $val eq "." ? 6 : $val =~ /[XHZU]/ ? 0 : $val =~ /[xhzu]/ ? 9 : ($SEQTEMP->[$i] ne "-" and $val ne "B") ? 2 : $val;
				$CT++ if $REFTEMP->[$i] eq "C" and $SEQTEMP->[$i] eq "T" and $BADTEMP->[$i] eq " ";
				$GA++ if $REFTEMP->[$i] eq "G" and $SEQTEMP->[$i] eq "A" and $BADTEMP->[$i] eq " ";
				#$tempstrand = $CT > $GA ? 0 : 16;
				push(@{$ref3}, $REFTEMP->[$i]);
				push(@{$seq3}, $SEQTEMP->[$i]);
				push(@{$myc3}, $myc4->[$i]);
				push(@{$con3}, $con2->[$i]);
				push(@{$bad3}, $BADTEMP->[$i]);
				$val = $val !~ /^[0-9]$/ ? -1 : $val;
				$valcg = $valcg !~ /^[0-9]$/ ? -1 : $valcg;
				$VAL2->[$i] =$val;
				$VALCG2->[$i] =$valcg;
			}
		}
		#print "AFTER: CT=$CT, GA=$GA\n";
	}
	my $output = $tempstrand == 0 ? $out{$chr}{OUTPOS} : $tempstrand == 16 ? $out{$chr}{OUTNEG} : die "Cannot determine tempstrand ($tempstrand)\n";
	my $outputCG = $tempstrand == 0 ? $out{$chr}{OUTPOSCG} : $tempstrand == 16 ? $out{$chr}{OUTNEGCG} : die "Cannot determine tempstrand ($tempstrand)\n";
	open (my $out, ">>", $output) or die "Cannot write to $output: $!\n";
	open (my $outCG, ">>", $outputCG) or die "Cannot write to $outputCG: $!\n";
	#print $out $read . " CT=$CT,GA=$GA,$strand->$tempstrand" . " " .  join("", @{$myc4}) . "\n";
	print $out "\"$read,CT=$CT,GA=$GA,$strand,$tempstrand\"\t" . join("\t", @{$VAL}) . "\n" if $tempstrand eq $strand;
	print $out "\"$read,CT=$CT,GA=$GA,$strand,$tempstrand\"\t" . join("\t", @{$VAL2}) . "\n" if $tempstrand ne $strand;
	#print $outCG $read . " CT=$CT,GA=$GA,$strand->$tempstrand" . " " .  join("", @{$myc4}) . "\n";
	print $outCG "\"$read,CT=$CT,GA=$GA,$strand,$tempstrand\"\t" . join("\t", @{$VALCG}) . "\n" if $tempstrand eq $strand;
	print $outCG "\"$read,CT=$CT,GA=$GA,$strand,$tempstrand\"\t" . join("\t", @{$VALCG2}) . "\n" if $tempstrand ne $strand;
#	if ($tempstrand ne $strand) {
	#	print $read . "\tCT=$CT,GA=$GA,$strand->$tempstrand" . "\t" .  join("", @{$myc4}) . "\n";
	#	print $read . "\tCT=$CT,GA=$GA,$strand->$tempstrand" . "\t" .  join("", @{$VAL2}) . "\n";
	#	print $read . "\tCT=$CT,GA=$GA,$strand->$tempstrand" . "\t" .  join("", @{$myc4}) . "\n";
	#	print $read . "\tCT=$CT,GA=$GA,$strand->$tempstrand" . "\t" .  join("", @{$VALCG2}) . "\n";
	#	die;
#	}
	print $outdebug "\n$YW$read$N\t$strand\t$pos\n";
	print $outdebug "NOCG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$ref3}) . "\n";
	print $outdebug "NOCG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$seq3}) . "\n";
	print $outdebug "NOCG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$myc4}) . "\n";
	print $outdebug "NOCG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$bad3}) . "\n";
	print $outdebug "NOCG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$VAL}) . "\n";
	print $outdebug "NOCG\t" . "CT=$CT, GA=$GA\n";
	print $outdebug "CG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$ref3}) . "\n";
	print $outdebug "CG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$seq3}) . "\n";
	print $outdebug "CG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$myc4}) . "\n";
	print $outdebug "CG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$bad3}) . "\n";
	print $outdebug "CG\t" . $read . "\t$strand->$tempstrand" . "\t" .  join("", @{$VALCG}) . "\n";
	print $outdebug "CG\t" . "CT=$CT, GA=$GA\n";
	$strand{$strand}{'same'} ++ if $strand eq $tempstrand;
	$strand{$strand}{'diff'} ++ if $strand ne $tempstrand;
	my $div = $CT+$GA == 0 ? 0 : int($CT/($CT+$GA)*1000+0.5)/10;
	my $typez = $strand eq $tempstrand ? "same" : "diff";
	push(@{$strand{$strand}{CT}{$typez}}, $div);
	push(@{$strand{$strand}{tot}{$typez}}, $CT+$GA);
	for (my $i = 0; $i < @{$ref2}; $i++) {
#	for (my $i = $seqborder0; $i < $seqborder1; $i++) {
		next if $i < $seqborder0 or $i >= $seqborder1;
		$count{Cref} ++ if $ref2->[$i] eq "C";
		$count{Gref} ++ if $ref2->[$i] eq "G";
		$count{Cseq} ++ if $seq2->[$i] eq "C";
		$count{Gseq} ++ if $seq2->[$i] eq "G";
		my $dinuc = $ref2->[$i] . $seq2->[$i];
		if ($con2->[$i] =~ /[hzxu]/) {
			$con2->[$i] = $LGN . $con2->[$i] . $N if not defined $bad2->[$i];
		}
		if ($dinuc =~ /[CG]/ and $ref2->[$i] eq $seq2->[$i]) {
			$count{$dinuc} ++;
			next;
		}
		elsif (defined $bad{$i} and $dinuc =~ /[CG]/) {
			$count{$dinuc . "bad"} ++;
			#next;
		}
		elsif ($dinuc eq "CT" or $dinuc eq "GA") {
			$count{$dinuc} ++;
			$ref2->[$i] = $LRD . $ref2->[$i] . $N;
			$seq2->[$i] = $LGN . $seq2->[$i] . $N;
		}
		elsif ($dinuc eq "TC" or $dinuc eq "AG") {
			$count{$dinuc} ++;
			$seq2->[$i] = $LRD . $seq2->[$i] . $N;
			$ref2->[$i] = $LGN . $ref2->[$i] . $N;
		}
		elsif ($ref2->[$i] =~ /^(C|G)$/) {
			$count{$dinuc} ++;
			$ref2->[$i] = $LRD . $ref2->[$i] . $N;
			$seq2->[$i] = $LPR . $seq2->[$i] . $N;
		}
		elsif ($seq2->[$i] =~ /^(C|G)$/) {
			$count{$dinuc} ++;
			$seq2->[$i] = $LRD . $seq2->[$i] . $N;
			$ref2->[$i] = $LPR . $ref2->[$i] . $N;
		}
	}
#	print join("", @{$VAL}) . "\n";
	#print join("", @{$VAL}) . "\n" if $strand eq 0;
	#print join("", @{$VAL}) . "\n" if $strand eq 16;
	($ref3, $seq3, $con3, $myc3, $bad3) = ();
	for (my $i = 0; $i < @{$ref2}; $i++) {
		$ref3->[$i] = $ref2->[$i] ne "-" ? $ref2->[$i] : " ";
		$seq3->[$i] = $ref2->[$i] ne "-" ? $seq2->[$i] : " ";
		$con3->[$i] = $ref2->[$i] ne "-" ? $con2->[$i] : " ";
		$myc3->[$i] = $ref2->[$i] ne "-" ? $myc2->[$i] : " ";
		$bad3->[$i] = $ref2->[$i] ne "-" ? $bad2->[$i] : " ";
	}
	if ($linecount >= 0) {
		print $outdebug ">$read\t$strand\t$pos\n";
		print $outdebug "BEFORE\n";
		print $outdebug join("", @{$ref2}) . "\n";	
		print $outdebug join("", @{$seq2}) . "\n";	
		print $outdebug join("", @{$con2}) . "\n";	
		print $outdebug join("", @{$myc2}) . "\n";	
		print $outdebug join("", @{$bad2}) . "\n";	
		print $outdebug "AFTER\n";
		print $outdebug join("", @{$ref3}) . "\n";	
		print $outdebug join("", @{$seq3}) . "\n";	
		print $outdebug join("", @{$con3}) . "\n";	
		print $outdebug join("", @{$myc3}) . "\n";	
		print $outdebug join("", @{$bad3}) . "\n";	
		foreach my $dinuc (sort keys %count) {
			next if $dinuc =~ /bad/ and $dinuc !~ /\-/; next if $dinuc =~ /(seq|ref)/;
			my ($nuc1, $nuc2) = split("", $dinuc);
			my $tot1 = $nuc1 eq "G" ? $count{Gref} : $nuc1 eq "C" ? $count{Cref} : 0;
			my $tot2 = $nuc2 eq "G" ? $count{Gseq} : $nuc2 eq "C" ? $count{Cseq} : 0;
			my $bad = defined $count{$dinuc . "bad"} ? $count{$dinuc . "bad"} : 0;
			print $outdebug "$dinuc: $LGN$count{$dinuc}$N+$LRD$bad$N (total ref=$LCY$tot1$N, total seq=$LPR$tot2$N)\n";
#			print $out "$dinuc: $LGN$count{$dinuc}$N+$LRD$bad$N (total ref=$LCY$tot1$N, total seq=$LPR$tot2$N)\n";# if $dinuc =~ /(CT|GA)/;
#			print $outCG "$dinuc: $LGN$count{$dinuc}$N+$LRD$bad$N (total ref=$LCY$tot1$N, total seq=$LPR$tot2$N)\n";# if $dinuc =~ /(CT|GA)/;
		}
	}
	#print $out "\n"; 
	close $out;
	#print $outCG "\n"; 
	close $outCG;
	($ref3, $seq3, $con3, $myc3, $bad3) = ();
	my $cons3 = join("", @{$con2}); my ($beg00, $end00) = $cons3 =~ /^([\-]+)([A-Z\.]+.+[A-Z\.])[\-]+$/;
	#print "BEG=$beg00\nEND=$end00\n";
	$beg00 = defined $beg00 ? length($beg00) : 0;
	$end00 = defined $end00 ? length($end00) + $beg00 : 0;
	$cons3 = "";
	my $cons4;
#	print "diff\n" if $strand ne $tempstrand;
	for (my $i = 0; $i < @{$ref2}; $i++) {
		if ($ref2->[$i] ne "-") {
			push(@{$ref3}, $ref2->[$i]);
			push(@{$seq3}, $seq2->[$i]);
			if ($tempstrand ne $strand) {
				$cons3 = $i < $beg00 ? " " : $i > ($end00-1) ? " " : $myc4->[$i] =~ /^[ \-]$/ ? "." : $myc4->[$i];
			}
			else {
				$cons3 = $con2->[$i];
			}
			push(@{$con3}, $cons3);
			if ($tempstrand eq $strand) {$cons3 = $myc2->[$i]}
			push(@{$myc3}, $cons3);
			push(@{$bad3}, $bad2->[$i]);
		}
	}
		print $outdebug "AFTER2\n";
		print $outdebug join("", @{$ref3}) . "\n";	
		print $outdebug join("", @{$seq3}) . "\n";	
		print $outdebug join("", @{$con3}) . "\n";	
		print $outdebug join("", @{$myc3}) . "\n";	
		print $outdebug join("", @{$bad3}) . "\n";	
	print $outdebug "strand=$strand, new strand=$tempstrand\n";# if $strand ne $tempstrand;
	$cons3 = join("", @{$con3}); $cons3 =~ s/[ ]+//g;
	$cons3 =~ s/^[\-]+//g;
	$cons3 =~ s/[\-]+$//g;
	$cons3 =~ s/[\-]/./g;
	#my ($read, $tempstrand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted) = split("\t", $line);
	$others = $tempstrand eq 0 ? "XR:Z:CT\tXG:Z:CT" : "XR:Z:CT\tXG:Z:GA";
#	print $outsam "$read\t$tempstrand\t$chr\t$pos\t$mapq\t$cigar\t$junk1\t$junk2\t$junk3\t$seqs\t$qual\t$junk4\t$junk5\tXM:Z:$cons3$others\n";
	print $outsam "$read\t$tempstrand\t$cons3\n";
#	last if $strand ne $tempstrand;
#	last if $linecount > 200;
}
close $outsam;
foreach my $strand (sort keys %strand) {
	my @types = ("same","diff");
	print "$strand: ";
	foreach my $type (@types[0..1]) {
		my $total = $strand{$strand}{$type}; $total = 0 if not defined $total;
		print "$type=$total,";
		my $CT = $strand{$strand}{CT}{$type};
		my ($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $CT) {
			$tmm = int(1000*mitochy::tmm(@{$CT})+0.5)/1000;
			$mean = int(1000*mitochy::mean(@{$CT})+0.5)/1000;
			$tmmse = int(1000*mitochy::tmmse(@{$CT})+0.5)/1000;
			$meanse = int(1000*mitochy::se(@{$CT})+0.5)/1000;
		}
		print "CT=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse, ";
		my $tot = $strand{$strand}{tot}{$type}; 
		($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $tot) {
			$tmm = int(1000*mitochy::tmm(@{$tot})+0.5)/1000;
			$mean = int(1000*mitochy::mean(@{$tot})+0.5)/1000;
			$tmmse = int(1000*mitochy::tmmse(@{$tot})+0.5)/1000;
			$meanse = int(1000*mitochy::se(@{$tot})+0.5)/1000;
		}
		print "tot=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse\n";
	}
}
exit 0;
if (defined $opt_C) {
	my @cons; my @tots;
	my @col = ($LCY, $LGN, $N, $YW, $LRD);
	my %cons2; my %tots2;
	for (my $i = 0; $i < @{$REF}; $i++) {
		my ($nuccol, $totcol);
		if (not defined $cons->[$i]) {
			push(@cons, " "); push(@tots, " "); next;
		}
		if (defined $cons->[$i] and defined $cons->[$i]{nuc2}) {
			my ($total2, $totcol2, $nuccol2); $total2 = 0;
			die if not defined $cons->[$i]{nuc2};
			my ($temp, $len);
			foreach my $len2 (sort {$cons->[$i]{nuc2}{$b}{tot} <=> $cons->[$i]{nuc2}{$a}{tot}} keys %{$cons->[$i]{nuc2}}) {
				$total2 = $cons->[$i]{nuc2}{$len2}{tot};
				$len = $len2;
				foreach my $nuc (sort {$cons->[$i]{nuc2}{$len2}{nuc}{$b} <=> $cons->[$i]{nuc2}{$len2}{nuc}{$a}} keys %{$cons->[$i]{nuc2}{$len2}{nuc}}) {
					my @nuc = split("", $nuc);
					for (my $j = 0; $j < $len; $j++) {
						$temp->[$j]{nuc}{$nuc[$j]} ++;
						$temp->[$j]{tot} ++;
					}
				}
				last;
			}
			my ($totcol3, $nuccol3);
			for (my $j = 0; $j < $len; $j++) {
				my $total = $temp->[$j]{tot}; $total = 0 if not defined $total;
				my $totcol2 = int($total/($linecount*0.2)); $totcol2 = @col-1 if $totcol2 >= @col; $totcol2 = $col[$totcol2] . $totcol2 . $N;
				foreach my $nuc (sort {$temp->[$j]{nuc}{$b} <=> $temp->[$j]{nuc}{$a}} keys %{$temp->[$j]{nuc}}) {
					my $nuctot = $temp->[$j]{nuc}{$nuc};
					$nuccol2 = int($nuctot/($total*0.2)); $nuccol2 = @col-1 if $nuccol2 >= @col; $nuccol2 = $col[$nuccol2] . $nuc . $N;
					last;
				}
				$totcol2 = " " if not defined $totcol2;
				$nuccol2 = " " if not defined $nuccol2;
				$totcol3 .= $totcol2;
				$nuccol3 .= $nuccol2;
			}
			$cons2{$i}{total} = $total2;
			$cons2{$i}{tots} = $totcol3;
			$cons2{$i}{cons} = $nuccol3;
		}
		if (defined $cons->[$i] and not defined $cons->[$i]{tot}) {
			push(@cons, " "); push(@tots, " "); next;
		}
		if (defined $cons->[$i] and not defined $cons->[$i]{nuc}) {
			push(@cons, " "); push(@tots, " "); next;
		}

		my $total = $cons->[$i]{tot}; $total = 0 if not defined $total;
		$totcol = int($total/($linecount*0.2)); $totcol = @col-1 if $totcol >= @col; $totcol = $col[$totcol] . $totcol . $N;
		foreach my $nuc (sort {$cons->[$i]{nuc}{$b} <=> $cons->[$i]{nuc}{$a}} keys %{$cons->[$i]{nuc}}) {
			my $nuctot = $cons->[$i]{nuc}{$nuc};
			$nuccol = int($nuctot/($total*0.2)); $nuccol = @col-1 if $nuccol >= @col; $nuccol = $col[$nuccol] . $nuc . $N;
			last;
		}
		$nuccol = " " if not defined $nuccol;
		$totcol = " " if not defined $totcol;
		push(@cons, $nuccol);
		push(@tots, $totcol);
	}

	my @refreal;
	for (my $i = 0; $i < @{$REF}; $i++) {
		push(@refreal, $REF->[$i]) if $REF->[$i] ne "-";
		my $curr = $i/10;
		print $curr . join("", (" ")x(10-length($curr))) if $i % 10 == 0;
	}
	print "\n";
	for (my $i = 0; $i < @{$REF}; $i++) {
		print $i%10;
	}
	print "\n";
	print join("", @refreal) . "\n";	
	print join("", @cons) . "\n";	
	print join("", @tots) . "\n";	

	foreach my $i (sort {$a <=> $b} keys %cons2) {
		my $total = $cons->[$i]{tot}; $total = 0 if not defined $total;
		print "$i: total=$cons2{$i}{total} (> $total) / $linecount\n$cons2{$i}{cons}\n$cons2{$i}{tots}\n\n";# if $total < $cons2{$i}{total};
	}
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

sub make_consensus {
	my ($ref2, $seq2, $seqborder0, $seqborder1, $cons) = @_;
	my $pos = $seqborder0; my $curr;
	for (my $i = $seqborder0; $i < $seqborder1; $i++) {
		my $ref = $ref2->[$i];
		my $seq = $seq2->[$i];
		if ($ref ne "-") {
			$cons->[$pos]{nuc}{$seq} ++;
			$cons->[$pos]{tot} ++;
			$pos ++;
			next;
		}
		my $temp; my $lastnuc;
		for (my $j = $i; $j < $seqborder1; $j++) {
			$temp .= $ref2->[$j];
			if ($ref2->[$j] ne "-") {
				my $len = length($temp);
				$lastnuc = $ref2->[$j];
				$cons->[$pos]{nuc2}{$len}{nuc}{$temp} ++;
				$cons->[$pos]{nuc2}{$len}{tot} ++;
				$cons->[$pos]{nuc}{$lastnuc} ++;
				$cons->[$pos]{tot} ++;
				$pos ++;
				$i = $j;
				last;
			}
		}
	}
	return($cons);
}

sub check_file {
	my ($file, $type) = @_;
	die "$type file $file does not exist!\n" if ex($file) == 0;
	die "$type file $file is empty!\n" if -s $file == 0;
	my $filetype = `file -b --mime-type $file`; chomp($filetype);
	my @line = ($file =~ /.bam$/ or $filetype !~ /text/) ? `samtools view $file | head -n 200` : `head -n 200 $file`;
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
sub parse_logFile {
	my ($logFile) = @_;
	my ($samFile, $seqFile, $genez, $threshold);
	if (-e $logFile) {
		my @line = `cat $logFile`;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i]; chomp($line);
			($samFile) = $line =~ /^.+SAM.+[ ]+:(.+\\.[bs]am)($|[ ]?.+\(will be generated\).+$)/i if $line =~ /^\d+\.[ ]+\-S .+SAM.+[ ]+:/;
			if ($line =~ /\- Running.+bismark_genome_preparation.+\-\-bowtie2 /) {
				($seqFile) = $line =~ /\- Running.+bismark_genome_preparation.+\-\-bowtie2 (.+)$/;
				$seqFile .= "/geneIndexes.fa";
			}
#			print "$line\n" if $i > 20 and $i < 50;
			if ($line =~ /2. Parsing in sequence for genes from sequence file/) {
				for (my $j = $i+1; $j < @line; $j++) {
					#print "\t$line\n";
					last if $line[$j] =~ /SUCCESS.+Sequence has been parsed from fasta file/;
					my ($gene, $length) = $line[$j] =~ /^.+genez=(.+) \(.+\) Length=(\d+)$/;
					die if not defined $gene or not defined $length;
					$genez->{$gene} = $length;
					print "parse_logFile: parsed gene=$LCY$gene$N=\n";
				}
			}
			($threshold) = $line =~ /^\d+\.[ ]+\-t .+Threshold.+[ ]:(\-?\d+\.?\d*)$/ if $line =~ /^\d+\.[ ]+\-t .+Threshold.+[ ]:(\-?\d+\.?\d*)$/;
			last if ex([$samFile, $seqFile]) == 1 and defined $threshold and defined $genez;
			last if $line =~ /SUCCESS.+Sequence has been parsed from fasta file/;
		}
	}
	if (not -e $logFile or not defined $samFile) {
		my @sam = <$folder/*.[sb]am>;
		die "Cannot find a .sam..bam file in folder $CY$folder$N\n\n" if @sam == 0;
		die "More than one .sam/.bam file in folder $CY$folder$N\n-" . join("\n-", @sam) . "\n\n" if @sam > 1 and not -e $logFile;
		$samFile = $sam[0] if not -e $logFile;
		$samFile = $sam[0] if not defined $samFile;
	}
	die "Please run footLoop.pl first before running this!\n" if not -e $logFile;
	die "Cant find threshold from $logFile!\n" if not defined $threshold;
	return($samFile, $seqFile, $genez, $threshold);
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
			$reftemp .= $ref2->[$j];
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
	#print join("", @bad) . "\n\n";
	#print $prints . "\n";
	return(\%bad);
}
sub parse {
	my ($line, $refs) = @_;
	my @refs = @{$refs};
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted) = split("\t", $line);
	my $con = getConv($converted);
	my @con;
	my %pos;
	my @seq = split("",$seqs);
	my ($num, $alp, $lengthseq) = mitochy::parse_cigar($cigar); die if not defined $num;
	my @num  = @{$num}; my @alp = @{$alp};
	my @ref0 = @refs[0..$pos-2];
	my @ref = @refs[$pos-1..@refs-1];
	my ($seq, $ref, $seqpos, $refpos) = (\@seq, \@ref, 0, 0);
	my $lengthref = @ref; my $insref = 0; my $insseq = 0;
	for (my $i = 0; $i < @num; $i++) {
		for (my $j = 0; $j < $num[$i]; $j++) {
			($ref) = $alp[$i] eq "I" ? ins($ref, $refpos, "-", "ref") : $ref;
			($seq) = $alp[$i] eq "D" ? ins($seq, $seqpos, "-", "seq") : $seq; 
			($con) = $alp[$i] eq "D" ? ins($con, $seqpos, "-", "con") : $con; 
			$refpos ++;
			$seqpos ++;
			$insref ++ if $alp[$i] eq "I";
			$insseq ++ if $alp[$i] eq "D";
			$pos{seq}{$seqpos-$insseq+$pos-1} = $refpos-1+$pos-1;
			$pos{ref}{$refpos-$insref+$pos-1} = $seqpos-1+$pos-1;
		}
	}
	my $refend = $refpos - $insref;
	print "$read\t$strand\t$pos\n";
#	print join("", (@ref0, @{$ref})) . "\n";
#	print join("", (("-") x ($pos-1))) . join("", @{$seq}) . join("", (("-")x($lengthref-$refend))) . "\n";
#	print join("", (("-") x ($pos-1))) . join("", @{$con}) . join("", (("-")x($lengthref-$refend))) . "\n";
	@ref = (@ref0, @{$ref});
	my $seqborder0 = $pos-1;
	my $seqborder1 = $pos - 1 + @{$seq};
	@seq = ((("-") x ($pos-1)), @{$seq}, (("-")x($lengthref-$refend)));
	@con = ((("-") x ($pos-1)), @{$con}, (("-")x($lengthref-$refend)));
	my ($myc, $printbug) = det_C_type(\@ref, \@seq, \@con, $strand, $seqborder0, $seqborder1);
#	print join("", @{$myc}) . "\n";
#	print "\n";
	#die if $strand eq 0;
#	print "$printbug\n";
	return(\@ref, \@seq, \@con, $myc, \%pos, $seqborder0, $seqborder1);
}
#-X CHG -H CHH -Z CG upper = non conv
sub det_C_type {
	my ($ref, $seq, $con, $strand, $seqborder0, $seqborder1) = @_;
	my @data;
	my $printbug = "";
	my ($nuc1, $nuc2, $convert) = $strand eq 0 ? ("C","G","T") : ("G","C","A");
	for (my $i = 0; $i < @{$ref}; $i++) {
		$printbug .= join("", (" ")x($i)) . $ref->[$i] if $i >= 230 and $i <= 280;
		if ($ref->[$i] ne $nuc1) {
			$data[$i] = " ";
		}
		else {
			if (not defined $data[$i]) {
				$data[$i] = "H";
				$printbug .= " (0: H)" if $i >= 230 and $i <= 280;
			}
		}
#		$data[$i] = "H" if $ref->[$i] eq $nuc1 and not defined $data[$i];
		my ($len1, $len2, $next1, $next2) = (@{$ref}-2, @{$ref}-3);
		for (my $j = $i+1; $j < @{$ref}-1; $j++) {
			$len1 = $j; 
			$next1 = $ref->[$j] if $ref->[$j] ne "-";
			$printbug .= "." if $ref->[$j] eq "-" and $i >= 230 and $i <= 280;
			$printbug .= "$ref->[$j]" if $ref->[$j] ne "-" and $i >= 230 and $i <= 280;
			last if defined $next1;
		}
		for (my $j = $len1+1; $j < @{$ref}-2; $j++) {
			$len2 = $j; 
			$next2 = $ref->[$j] if $ref->[$j] ne "-";
			$printbug .= "." if $ref->[$j] eq "-" and $i >= 230 and $i <= 280;
			$printbug .= "$LGN$ref->[$j]$N" if $ref->[$j] ne "-" and $i >= 230 and $i <= 280;
			last if defined $next2;
		}
		my ($pos1, $pos2) = $strand eq 0 ? ($i,$i) : ($len1,$len2); #pos1 = C>G<, #pos2 = CH>G<
		if ($i < @{$ref}-1) {
			my $dinuc = $ref->[$i] . $ref->[$len1];
			$data[$pos1] = "Z" if $dinuc eq "CG"; 
			$printbug .= " ($pos1: Z)" if $dinuc eq "CG" and $i >= 230 and $i <= 280;
			$printbug .= " ($pos1: NOTCG)" if $dinuc ne "CG" and $i >= 230 and $i <= 280;
		}
		if ($i < @{$ref}-2) {
			my $dinuc = $ref->[$i] . $ref->[$len1];
			my $trinuc = $dinuc . $ref->[$len2];
			$printbug .= ">$data[$pos2]<" if defined $data[$pos2] and $i >= 230 and $i <= 280;
			$printbug .= ">UNDEF $trinuc<" if not defined $data[$pos2] and $i >= 230 and $i <= 280;
#			if (((defined $data[$pos2] and $data[$pos2] eq "H") or 
			if (not defined $data[$pos2] and $trinuc =~ /C[GACTNU]G/) {
				$data[$pos2] = "X";
			}
			elsif ($dinuc ne "CG" and $trinuc =~ /C[ACTNUG]G/ and (not defined $data[$pos2] or (defined $data[$pos2] and $data[$pos2] ne "Z"))) {
				$data[$pos2] = "X";
				$printbug .= "$pos2: (X)" if $i >= 230 and $i <= 280;
			}
			$printbug .= " ($pos2: NOTCHG)" if $i >= 230 and $i <= 280 and (not defined $data[$pos2] or $data[$pos2] ne "X");
		}
		$printbug .= "\n" if $i >= 230 and $i <= 280;
	}
	for (my $i = 0; $i < @{$ref}; $i++) {
		if ($i < $seqborder0 or $i >= $seqborder1) {$data[$i] = " "; next}
		if ($ref->[$i] . $seq->[$i] eq $nuc1 . $convert) {
			$data[$i] = lc($data[$i]);
		}
		elsif ($ref->[$i] . $seq->[$i] eq $nuc1 . $nuc1) {
			$data[$i] = ($data[$i]);
		}
		else {
			$data[$i] = ($ref->[$i] eq $nuc1 and $seq->[$i] ne $convert and $data[$i] eq "Z") ? "c" : $ref->[$i] eq $nuc1 ? "n" : $data[$i];
		}
		if ($data[$i] ne $con->[$i]) {
			#$data[$i] = $LRD . $data[$i] . $N;
		}
		#print "$i-2 to $i+2: " . join("", @{$ref}[$i-2..$i+2]) . "\n" if not defined $data[$i];
		die if not defined $data[$i];
	}
	$printbug = "";
	return(\@data, $printbug);
}
sub ins {
	my ($arr, $pos, $ins, $type) = @_;
	my @arr = @{$arr};
	my @arr0 = @arr[0..$pos-1];
	my @arr1 = @arr[$pos..@arr-1];
	#print "$type $pos: " . join("", @arr0) . "\n";
	die if @arr0 == 0;
	@arr = (@arr0, $ins, @arr1);
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
__END__
	my ($converted) = $arr[13] =~ /^.+:([\.A-Z]+)$/i; my $length2 = length($converted);
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
	
	my %count;
	my $ref2; my $seq2; my ($begseq, $endseq) = (0,0); my ($begref, $endref) = ($pos-1,$pos-1);
#	my $refprint = substr($ref, $pos-1, length($ref)-$pos); $refprint =~ s/A/${LPR}A$N/g; $refprint =~ s/C/${CY}C$N/g; $refprint =~ s/G/${LGN}G$N/g;$refprint =~ s/T/${YW}T$N/g;
	my $refprint = $ref; $refprint =~ s/A/${LPR}A$N/g; $refprint =~ s/C/${CY}C$N/g; $refprint =~ s/G/${LGN}G$N/g;$refprint =~ s/T/${YW}T$N/g;
	my $seqprint = join("", ("-") x ($pos-1)) . $seq; $seqprint =~ s/A/${LPR}A$N/g; $seqprint =~ s/C/${CY}C$N/g; $seqprint =~ s/G/${LGN}G$N/g;$seqprint =~ s/T/${YW}T$N/g;
	my $cigprint = join("", ("-") x ($pos-1));
	my @seq2print;# = join("", ("-")x($pos-1));
	my @ref2print;# = split("", substr($ref,0, $pos-1));
	for (my $i = 0; $i < 10; $i++) {
		$count{mis}{ID}{$i} = 1;
	}
	for (my $i = length($seq)-20; $i < length($seq); $i++) {
		$count{mis}{ID}{$i} = 1;
	}

	#print "$refprint\n$seqprint\n\n";
	while ($cigar =~ /\d+[A-Z]/g) {
		my ($prev, $curr, $next) = ($`, $&, $');
		my ($num, $cig) = $curr =~ /^(\d+)([A-Za-z]+)$/;
		$count{$cig} += $num/$length * 100;
		($begref, $begseq) = ($endref, $endseq);
		$endref = $cig =~ /^[MD]$/ ? $endref + $num : $endref;
		$endseq = $cig =~ /^[MI]$/ ? $endseq + $num : $endseq;
		my $refadd = ($cig eq "I" or $begref > length($ref)) ? join("", ("-") x ($num)) : join("", @ref[$begref..$endref-1]);
		my $seqadd = ($cig eq "D" or $begseq > length($seq)) ? join("", ("-") x ($num)) : join("", @seq[$begseq..$endseq-1]);
#		my $cigadd = ($cig eq "I") ? join("", ("I") x ($num)) : $cig eq "D" ? join("", ("D") x ($num)) : "";
		my $cigadd = ($cig eq "I") ? join("", (" ") x ($num)) : $cig eq "D" ? join("", (" ") x ($num)) : "";
		my $ref2printadd = ($cig eq "I") ? join("", ("-") x ($num)) : $cig =~ /^[D]$/ ? join("", @ref[$begref..$endref-1]) : join("", @ref[$begref..$endref-1]);
		push(@ref2print, split("", $ref2printadd));
		my $seq2printadd = ($cig =~ /^[I]$/) ? join("", @seq[$begseq..$endseq-1]) : $cig eq "D" ? join("", ("-") x ($num)) : join("", @seq[$begseq..$endseq-1]);
		push(@seq2print, split("", $seq2printadd));
		if ($cig =~ /^[ID]$/) {
			for (my $i = 0; $i < length($refadd); $i++) {
				$count{mis}{ID}{$begseq+$i} = 1 if $cig eq "I" or $cig eq "D";
			}
		}
		if ($cig eq "M") {
			for (my $i = 0; $i < length($refadd); $i ++) {
				my $ref3 = substr($refadd, $i, 1);
				my $seq3 = substr($seqadd, $i, 1);
				$count{m} += 1/$length * 100;
				$cigadd .= $ref3 eq $seq3 ? " " : "m";
				#push(@ref2print, $ref3);
				#push(@seq2print, $seq3);
				#$ref2printadd = ($ref3 !~ /^[CG]$/ and $seq3 !~ /^[CG]$/) ? $ref3 : $ref3 eq $seq3 ? "$ref3" : $ref3 =~ /^[CG]$/ ? "$LRD$ref3$N" : "$LGN$ref3$N";
				#push(@{@ref2print, $ref2printadd);
				#$seq2printadd = ($ref3 !~ /^[CG]$/ and $seq3 !~ /^[CG]$/) ? $seq3 : $ref3 eq $seq3 ? "$seq3" : $seq3 =~ /^[CG]$/ ? "$LRD$seq3$N" : "$LGN$seq3$N";
				#push(@{@seq2print, $seq2printadd);
				$count{mis}{CT}{$begseq+$i} = 1 if $ref3 eq "C" and $seq3 eq "T";
				$count{mis}{TC}{$begseq+$i} = 1 if $ref3 eq "T" and $seq3 eq "C";
				$count{mis}{GA}{$begseq+$i} = 1 if $ref3 eq "G" and $seq3 eq "A";
				$count{mis}{AG}{$begseq+$i} = 1 if $ref3 eq "A" and $seq3 eq "G";
				$count{mis}{GN}{$begseq+$i} = 1 if $ref3 eq "G" and $seq3 =~ /^[CT]$/;
				$count{mis}{NG}{$begseq+$i} = 1 if $seq3 eq "G" and $ref3 =~ /^[CT]$/;
				$count{mis}{CN}{$begseq+$i} = 1 if $ref3 eq "C" and $seq3 =~ /^[GA]$/;
				$count{mis}{CG}{$begseq+$i} = 1 if $seq3 eq "C" and $ref3 =~ /^[GA]$/;
			}	
		}
		$ref2 .= $refadd;
		$seq2 .= $seqadd;
		$cigprint .= $cigadd;
	}
	foreach my $cig (sort keys %count) {
		next if length($cig) != 1;
		$count{$cig} = int($count{$cig} * 10+0.5)/10;
		push(@{$data{$cig}}, $count{$cig});
	}
	push(@{$data{conv}}, $conv);
	push(@{$data{nonc}}, $nonc);
	push(@{$data{notc}}, $notc);
	push(@{$data{totl}}, $length);
	my $linecount = (@{$data{conv}});
	print "Done $linecount\n" if $linecount % 100 == 0;
	die "Died at count mis $count{mis}\n" if $count{mis} =~ /^\d+$/;
	for (my $pos = 0; $pos < @ref2print; $pos++) {
		#foreach my $pos (keys %{$count{mis}{$key}}) {
		my $ref3 = $ref2print[$pos];
		my $seq3 = $seq2print[$pos];
		my $checks = 0;
		foreach my $key (sort keys %{$count{mis}}) {
			next if $key eq "ID";
			next if not defined $count{mis}{$key}{$pos};
			for (my $k = -3; $k <= 3; $k++) {
				$checks = 1 if (defined $count{mis}{ID}{$pos+$k});# or defined $count{mis}{ID}{$pos-2} or defined $count{mis}{ID}{$pos+1} or defined $count{mis}{ID}{$pos+2}) {
				last if $checks eq 1;
			}
			if ($checks eq 1) {
				$count{mis}{$key}{$pos} -= 1;
=comment
				$count{mis}{CT}{$pos} -= 1 if $ref3 eq "C" and $seq3 eq "T";
				$count{mis}{TC}{$pos} -= 1 if $ref3 eq "T" and $seq3 eq "C";
				$count{mis}{GA}{$pos} -= 1 if $ref3 eq "G" and $seq3 eq "A";
				$count{mis}{AG}{$pos} -= 1 if $ref3 eq "A" and $seq3 eq "G";
				$count{mis}{GN}{$pos} -= 1 if $ref3 eq "G" and $seq3 =~ /^[CT]$/;
				$count{mis}{NG}{$pos} -= 1 if $seq3 eq "G" and $ref3 =~ /^[CT]$/;
				$count{mis}{CN}{$pos} -= 1 if $ref3 eq "C" and $seq3 =~ /^[GA]$/;
				$count{mis}{CG}{$pos} -= 1 if $seq3 eq "C" and $ref3 =~ /^[GA]$/;
=cut
				next;
			}
		}
		#($ref2print[$pos], $seq2print[$pos]) = ("${LRD}$ref3$N", "${LPR}$seq3$N") if "$ref3$seq3" =~ /^[CG].$/;
		#($ref2print[$pos], $seq2print[$pos]) = ("${LPR}$ref3$N", "${LRD}$seq3$N") if "$ref3$seq3" =~ /^.[CG]$/;
		($ref2print[$pos], $seq2print[$pos]) = ("${LRD}$ref3$N", "${LGN}$seq3$N") if "$ref3$seq3" =~ /^(CT|CA)$/;
		($ref2print[$pos], $seq2print[$pos]) = ("${LGN}$ref3$N", "${LRD}$seq3$N") if "$ref3$seq3" =~ /^(TC|AC)$/;
		($ref2print[$pos], $seq2print[$pos]) = ("${LRD}$ref3$N", "${LGN}$seq3$N") if "$ref3$seq3" =~ /^(GA|GT)$/;
		($ref2print[$pos], $seq2print[$pos]) = ("${LGN}$ref3$N", "${LRD}$seq3$N") if "$ref3$seq3" =~ /^(AG|TG)$/;
		#$ref2print[$pos] = ($ref3 !~ /^[CG]$/ and $seq3 !~ /^[CG]$/) ? $ref3 : $ref3 eq $seq3 ? "$ref3" : $ref3" =~ /^[CG]$/ ? "$LRD$ref3$N" : "$LGN$ref3$N";
		#$seq2print[$pos] = ($ref3 !~ /^[CG]$/ and $seq3 !~ /^[CG]$/) ? $seq3 : $ref3 eq $seq3 ? "$seq3" : $seq3" =~ /^[CG]$/ ? "$LRD$seq3$N" : "$LGN$seq3$N";
		#push(@{@ref2print, $ref2printadd);
	}
	foreach my $key (sort keys %{$count{mis}}) {
		next if $key eq "ID";
		foreach my $pos (sort {$a <=> $b} keys %{$count{mis}{$key}}) {
#			$count{tot}{$key} += (defined $count{mis}{ID}{$pos-1} or defined $count{mis}{ID}{$pos-2} or defined $count{mis}{ID}{$pos+1} or defined $count{mis}{ID}{$pos+2}) ? 0 : 1;
		}
	}
	my @cigs = qw(CT GA TC AG CN GN);
	foreach my $key (@cigs[0..@cigs-1]) {
		next if $key eq "ID";
		my $total = (keys %{$count{mis}{$key}}); $total = 0 if not defined $total;
		foreach my $pos (sort {$a <=> $b} keys %{$count{mis}{$key}}) {
			$count{tot}{$key} += (defined $count{mis}{ID}{$pos-1} or defined $count{mis}{ID}{$pos-2} or defined $count{mis}{ID}{$pos+1} or defined $count{mis}{ID}{$pos+2}) ? 0 : 1;
		}
		$count{tot}{$key} = 0 if not defined $count{tot}{$key};
		print "$key=$count{tot}{$key}/$total; ";
	}
	print "\n";
	my $seq2print = join("", @seq2print);
	my $ref2print = join("", @ref2print);
	$seq2print = join("", ("-")x($pos-1)) . join("", @seq2print);
	$ref2print = substr($ref,0, $pos-1) . join("", @ref2print);
	#my $ref2print = substr($ref, 0, $pos-1) . $ref2; $ref2print =~ s/C/${LRD}C$N/g; $ref2print =~ s/G/${LGN}G$N/g; #$ref2print =~ s/G/${LGN}G$N/g;$ref2print =~ s/T/${YW}T$N/g;
	#my $seq2print = join("", ("-")x($pos-1)) . $seq2; $seq2print =~ s/C/${LRD}C$N/g; $seq2print =~ s/G/${LGN}G$N/g;# $seq2print =~ s/G/${LGN}G$N/g;$seq2print =~ s/T/${YW}T$N/g;
	print "$strand\n$ref2print\n$seq2print\n$cigprint\n\n";die if $linecount > 10;
}
close $in1;

my ($conv) = mitochy::tmm(@{$data{conv}});
my ($nonc) = mitochy::tmm(@{$data{nonc}});
my ($notc) = mitochy::tmm(@{$data{notc}});
my ($totl) = mitochy::tmm(@{$data{totl}});
foreach my $cig (sort keys %data) {
	my $cigar2 = $cig eq "M" ? "Mat" : $cig eq "m" ? "Mis" : $cig eq "D" ? "Del" : $cig eq "I" ? "Ins" : "UNK";
	next if length($cig) > 1;
	my ($tmm) = mitochy::tmm(@{$data{$cig}});
	my $num = int($tmm / 100 * $totl * 10 +0.5)/10;
	$tmm = int($tmm * 10+0.5)/10;
	print "$cig\t$tmm ($num)\n";
}
print "total=$totl, conv=$conv, non-conv=$nonc, notc=$notc\n";


=comment
	my $prev = $pos;
	print join("", (" ")x($pos-1));
	for (my $i = $pos-1; $i < @ref1; $i++) {
		my $curr = $poz{ref}{$i}; $curr = $i-1 if not defined $curr;
#		print "i=$i, curr=$curr\n";
		print join("", (".")x($curr-1-$prev)) . $ref1[$i];
#		print join("", (" ")x($curr)) . $ref1[$i] . "\n";
		$prev = $curr;
	}
	print "\n\n";
=cut
