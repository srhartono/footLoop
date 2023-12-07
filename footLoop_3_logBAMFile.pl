#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite; use myFootLib;
use vars qw($opt_v $opt_b $opt_L $opt_g $opt_i $opt_o $opt_B $opt_q $opt_x $opt_y $opt_r $opt_O);
getopts("vb:L:g:i:o:B:x:y:r:O:");

my ($BAMFile, $minReadL, $seqFile, $geneIndexFile, $outDir, $filteredDir) = ($opt_b, $opt_L, $opt_g, $opt_i, $opt_o, $opt_O);
die "\nusage: $YW$0$N -r <readFile.fq.gz> -b $LCY<BAMFile>$N -L $LGN<minReadL>$N -g $CY<seqFile>$N -i $LCY<geneIndexFile>$N\n\n" unless defined $BAMFile and defined $minReadL and defined $seqFile and defined $geneIndexFile and -e $BAMFile and -e $seqFile and -e $geneIndexFile;
my ($readFilename)  = getFilename($opt_r, "full");
my $minMapQ    = (not defined($opt_q)) ? 0  : $opt_q;
($minReadL)   = $opt_L =~ /p$/i ? $opt_L =~ /^(.+)p$/i : $opt_L;
my $bufferL    = (not defined($opt_x)) ? 0  : $opt_x;
my $bufferR    = (not defined($opt_y)) ? 0  : $opt_y;
my $bufferLen  = $bufferR - $bufferL;

makedir($outDir) if not -d $outDir;

my ($folder1, $fileName1) = getFilename($BAMFile, "folderfull");
open (my $outLog, ">", "$outDir/$fileName1.log.txt") or die "Failed to write to $BAMFile.log.txt: $!\n";
my $SEQ = parse_fasta($seqFile, $outLog, $minReadL);
#my $filteredDir = $BAMFile . "_splitBAM";
parse_BAMFile($BAMFile, $seqFile, $filteredDir, $outLog);
die "good debug\n";

LOG($outLog, "outDir=$outDir\n");
close $outLog;

system("echo DONE > $outDir/$fileName1.logBAMFile.done");

sub parse_fasta {
   my ($seqFile, $outLog, $minReadL) = @_;
   #LOG($outLog, "\n\tc. Parsing in gene sequence and infos from seqFile=$CY$seqFile$N\n");
   open(my $SEQIN, "<", $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
   my $fasta = new FAlite($SEQIN);
   my $linecount = 0;
   #LOG($outLog, "\t\t${GR}From fasta file:$N\n");
   while (my $entry = $fasta->nextEntry()) {
      $linecount ++;
      my $def = $entry->def; $def =~ s/^>//;
      if ($def =~ /^>.+::.+:\d+\-\d+$/) { #newer version of fastaFromBed
         ($def) = $def =~ /^(.+)::.+:\d+\-\d+$/;
      }
      my $gene = uc($def);
      my $seqz = uc($entry->seq);
      $SEQ->{$gene}{seq}       = [split("", $seqz)];
      $SEQ->{$gene}{geneL}     = scalar(@{$SEQ->{$gene}{seq}}) - $bufferLen;
      $SEQ->{$gene}{geneL}     = scalar(@{$SEQ->{$gene}{seq}}) if $SEQ->{$gene}{geneL} <= 0;
      $SEQ->{$gene}{minReadL}  = $minReadL; #(defined $minReadL and $opt_L =~ /p$/i) ? int(0.5+$SEQ->{$gene}{geneL} * $minReadL / 100) : $minReadL;
      $SEQ->{$gene}{total}     = 0;
      $SEQ->{$gene}{badlength} = 0;
      $SEQ->{$gene}{lowq}      = 0;
      $SEQ->{$gene}{used}      = 0;
      $SEQ->{$gene}{pos}       = 0;
      $SEQ->{$gene}{neg}       = 0;
      $SEQ->{$gene}{orig}      = $def;
    #  LOG($outLog, "\t\t$GR$linecount:gene=$gene,length=$SEQ->{$gene}{geneL}$N\n");
   }
   close $SEQIN;
   LOG($outLog, date() . "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n");
   return($SEQ);
}

sub parse_BAMFile {
	my ($BAMFile, $seqFile, $filteredDir, $outLog) = @_;
	my ($BAMFileName) = getFilename($BAMFile);

	LOG($outLog, "\ta. Parsing BAM file $CY$BAMFile$N and getting only high quality reads\nfilename=$BAMFileName\n\n");
	open(my $notused, ">", "$outDir/.$readFilename.notused") or (LOG($outLog, "Cannot open $outDir/.$readFilename.notused: $!\n") and exit 1);
	my $BAM;
	#open($BAM, $BAMFile) or (LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $BAMFile: $!") and exit 1 if $BAMFile =~ /.sam$/);
	open($BAM, "samtools view $BAMFile|") or (LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $BAMFile: $!") and exit 1) if $BAMFile =~ /.bam$/;

	## Some stats
	my $linecount   = 0;
	my $BAMStats = makehash(['total','used','diffgene','lowq','badlength']);

	while(my $line = <$BAM>) {
		if ($linecount == 5) {last;}#DIELOG($outLog, "DEBUG GOOD\n");}
		$linecount ++;
		chomp($line);

		#####################
		# 1. Parse BAM line #
		#####################
		my @arr = split("\t", $line);
		LOG($outLog, "\t$YW$BAMFile$N: Done $linecount\n") if $linecount % 5000 == 0;

		# a. Total column must be 14 or skipped (e.g. BAM header)
		if (@arr < 14) {
			LOG($outLog, "$LGN$linecount$N: BAM header detected (#column < 14). LINE:\n\t$line\n");
			next;
		}

		my ($eval, $evalPrint) = myeval(\@arr);
		my ($read, $readStrand, $gene, $readPos, $mapQ, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $tags, $junk4, $readMethCall) = @arr;
		my $seqLen = length($seqs);
		$gene = uc($gene);
		check_BAM_field($BAMFile, $outLog, @arr);
		LOG($outLog, "\t$BAMFile: length of gene $gene undef\n") and exit 1 if not defined($SEQ->{$gene}{geneL});

		$BAMStats->{total} ++;
		$BAMStats->{gene}{$gene} ++;

		# Filter out reads which gene isn't in genomeFile and indexFile
		if (not defined $SEQ->{$gene}) {
			LOG($outLog, "$LGN$linecount$N: read $LCY$read$N is mapped to gene $LGN$gene$N is not in geneIndex! LINE:\n\t$line\n") if not defined $SEQ->{$gene};
			next;
		}

		$SEQ->{$gene}{total} ++;
		my ($readname) = length($read) >= 20 ? $read =~ /(.{20})$/ : $read; $readname = "..." . $readname  if length($read) >= 20; $readname = $read if not defined $readname;
		LOG($outLog, "\tExample at read $BAMStats->{total}: name=$CY$readname$N\tstrand=$CY$readStrand$N\tchr/gene=$CY$gene$N\tpos=$CY$readPos$N\tmapQ=$CY$mapQ$N\n\n") if $BAMStats->{total} == 1;
		LOG($outLog, "\tDone $GN$BAMStats->{total}$N\n") if $BAMStats->{total} % 2000 == 0;


		# a. Filter out reads with mapping quality=$mapQ less than threshold quality=$opt_q
		if($mapQ < $minMapQ) {
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'lowq');
			print $notused "\t$CY$read$N quality ($CY$mapQ$N) is less than $CY$opt_q$N\n";
			$BAMStats->{lowq} ++;
			$SEQ->{$gene}{lowq} ++;
			next;
		}

		# b. Filter out reads with read length=$seqLen
		elsif (parse_cigar($cigar, "len") < $SEQ->{$gene}{minReadL}) {
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'badlength');
			my $cigarLen = parse_cigar($cigar, "len");
			print $notused "\t$CY$read$N length of seq (cigarLen = $LGN$cigarLen$N, seqLen = $CY$seqLen$N) is less than $SEQ->{$gene}{minReadL} bp (length of original sequence is ($CY" . $SEQ->{$gene}{geneL} . "$N)!\n";
		}

		# c. Filter out reads with more than 5% indel
		elsif (count_indel($cigar) > 5) {
			($SEQ, $BAMStats) = filter_BAM_read($SEQ, $BAMStats, $readname, $gene, 'indel');
			my ($indelperc) = count_indel($cigar);

			print $notused "\t$CY$read$N has high indel ($indelperc)!\n";
		}

		# c. If read strand is 0 or 16 store it
		elsif ($readStrand == 0 || $readStrand == 16)   {
			$SEQ->{$gene}{read}{$read} = $readStrand;
		}

		# d. Otherwise put it into notused
		else {
			print $notused "\t$CY$read$N isn't used for some reason\n";
		}
	}
	LOG($outLog, "\n\tb.Logging BAM file $CY$BAMFile$N\n");
	log_BAMFile($SEQ, $BAMFile, $BAMStats, $filteredDir, $outLog);
}

sub check_BAM_field {
   my ($BAMFile, $outLog, @arr) = @_;
   my $check = 0;
   for (my $i = 0; $i < @arr; $i++) {
      if ( not defined $arr[$i]) {
         LOG($outLog, date() . "footLoop.pl::parse_BAMFile::check_BAM_field $BAMFile: column $i undefined\n");
         $check = 1;
      }
   }
   DIELOG($outLog, "footLoop.pl::parse_BAMFile::check_BAM_field: $BAMFile: Something is wrong with the bam file!\n") if $check != 0;
}
sub count_indel {
   my ($cigar) = @_;
   my ($nums, $alps, $lenz) = parse_cigar($cigar);
   my ($ins, $del, $mat, $oth, $tot) = (0,0,0,0,0);
   for (my $i = 0; $i < @{$nums}; $i++) {
      my $alp = $alps->[$i];
      my $num = $nums->[$i];
      $mat += $num if $alp eq "M";
      $ins += $num if $alp eq "I";
      $del += $num if $alp eq "D";
      $oth += $num if $alp !~ /^(M|I|D)$/;
      $tot += $num;
   }
   my $matperc = $tot == 0 ? 0 : int($mat/$tot*1000+0.5)/10;
   my $insperc = $tot == 0 ? 0 : int($ins/$tot*1000+0.5)/10;
   my $delperc = $tot == 0 ? 0 : int($del/$tot*1000+0.5)/10;
   my $othperc = $tot == 0 ? 0 : int($oth/$tot*1000+0.5)/10;
   return $insperc+$delperc;
}



sub log_BAMFile {
	my ($SEQ, $BAMFile, $BAMStats, $filteredDir, $outLog) = @_;
	foreach my $genez (sort keys %{$SEQ}) {
		my $outTXTFilePos  = "$filteredDir/$genez\_Pos.filtered.gz";
		my $outTXTFileNeg  = "$filteredDir/$genez\_Neg.filtered.gz";
		#my $outTXTFileUnk  = "$filteredDir/$genez\_Unk.filtered.gz";
		open ($SEQ->{$genez}{outTXTPos}, "| gzip -f > $outTXTFilePos");
		open ($SEQ->{$genez}{outTXTNeg}, "| gzip -f > $outTXTFileNeg");
		#open ($SEQ->{$genez}{outTXTUnk}, "| gzip -f > $outTXTFileUnk");
	}
	my ($BAMFileName) = getFilename($BAMFile);

	my $skipped = 0; my ($passedFilterP, $passedFilterN) = (0,0);
	open (my $inBAMFix, "zcat $filteredDir/$BAMFileName.fixed.gz|") or LOG($outLog, "Failed to open $filteredDir/$BAMFileName.fixed.gz: $!\n") and exit 1;
	while (my $line = <$inBAMFix>) {
		chomp($line);
		my ($read, $type, $oldStrand, $strand, $genez, $pos, $info) = split("\t", $line);
		if (not defined $SEQ->{$genez} or not defined $info) {
			LOG($outLog, "\tfootLoop.pl: ERROR in $LCY$filteredDir/$BAMFileName.fixed.gz$N: gene=$genez but \$SEQ->{\$genez} is not defined!line=\n$line\n\n") if not defined $SEQ->{$genez};
			DIELOG($outLog, "\tfootLoop.pl: ERROR in $LCY$filteredDir/$BAMFileName.fixed.gz$N: gene=$genez and seq genez is $SEQ->{$genez} but info is not defined!line=\n$line\n\n") if defined $SEQ->{$genez} and not defined $info;
			$skipped ++;
			next;
		}
		if (not defined $SEQ->{$genez}{read}{$read}) {
			next;
		}
		#my ($CT0, $CC0, $GA0, $GG0, $CT1, $CC1, $GA1, $GG1) = split(",", $info);
		#if ($type eq "6_BOTH" or $type eq "3_NONE") {
		#	print {$SEQ->{$genez}{outTXTUnk}} "$read\tBP\t$pos\n" if $strand eq 0;
		#	print {$SEQ->{$genez}{outTXTUnk}} "$read\tBN\t$pos\n" if $strand eq 16;
		#	$SEQ->{$genez}{unkpos} ++ if $strand == 0;
		#	$SEQ->{$genez}{unkneg} ++ if $strand == 16;
		#	$passedFilterP ++ if $strand == 0;
		#	$passedFilterN ++ if $strand == 16;
		#}
		#else {
			print {$SEQ->{$genez}{outTXTNeg}} "$read\tFN\t$pos\n" if $strand eq 16;
			print {$SEQ->{$genez}{outTXTPos}} "$read\tFP\t$pos\n" if $strand eq 0;
			$SEQ->{$genez}{pos} ++ if $strand == 0;
			$SEQ->{$genez}{neg} ++ if $strand == 16;
			$passedFilterP ++ if $strand == 0;
			$passedFilterN ++ if $strand == 16;
		#}
		$BAMStats->{used} ++;
		$SEQ->{$genez}{used} ++;
		$SEQ->{$genez}{posneg} ++ if $strand eq 16 and $oldStrand eq 0;
		$SEQ->{$genez}{negpos} ++ if $strand eq 0 and $oldStrand eq 16;
	}
	close $inBAMFix;
	#LOG($outLog, "\t$LCY$filteredDir/$BAMFileName.fixed$N: skipped = $LGN$skipped$N\n");

	#LOG($outReadLog, "footLoop.pl,read_passed_filter,header\ttotal\tpositive\tnegative\n");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter,record\t$BAMStats->{total}\t$passedFilterP\t$passedFilterN\n");

	LOG($outLog, "
Reads that passed filters:
Positive: $passedFilterP
Negative: $passedFilterN
Total   : $BAMStats->{total};

Per Gene:
	","NA");

	my $zero = "";
	foreach my $gene (sort keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		foreach my $key (@key) {
			$SEQ->{$gene}{$key} = 0 if not defined $SEQ->{$gene}{$key};
		}
	}
	foreach my $gene (sort {$SEQ->{$b}{total} <=> $SEQ->{$a}{total}} keys %{$SEQ}) {
		my @key = qw(posneg negpos unkpos unkneg pos neg used total badlength lowq orig);
		my $outTXTFilePos  = "$filteredDir/$gene\_Pos.filtered.gz"; #system("/bin/rm $outTXTFilePos") if -e $outTXTFilePos and -s $outTXTFilePos == 0;
		my $outTXTFileNeg  = "$filteredDir/$gene\_Neg.filtered.gz"; #system("/bin/rm $outTXTFileNeg") if -e $outTXTFileNeg and -s $outTXTFileNeg == 0;
		#my $outTXTFileUnk  = "$filteredDir/$gene\_Unk.filtered"; system("/bin/rm $outTXTFileUnk") if -e $outTXTFileUnk and -s $outTXTFileUnk == 0;

		$zero .= "$gene ($SEQ->{$gene}{total})\n" and next if $SEQ->{$gene}{total} <= 10;
		$SEQ->{$gene}{total} = 0 if not defined $SEQ->{$gene}{total};
		$SEQ->{$gene}{orig} = 0 if not defined $SEQ->{$gene}{orig};
		$SEQ->{$gene}{negpos} = 0 if not defined $SEQ->{$gene}{negpos};
		$SEQ->{$gene}{posneg} = 0 if not defined $SEQ->{$gene}{posneg};
		$SEQ->{$gene}{unkpos} = 0 if not defined $SEQ->{$gene}{unkpos};
		$SEQ->{$gene}{unkneg} = 0 if not defined $SEQ->{$gene}{unkneg};
		$SEQ->{$gene}{used} = 0 if not defined $SEQ->{$gene}{used};
		$SEQ->{$gene}{total} = 0 if not defined $SEQ->{$gene}{total};
		$SEQ->{$gene}{badlength} = 0 if not defined $SEQ->{$gene}{badlength};
		$SEQ->{$gene}{indel} = 0 if not defined $SEQ->{$gene}{indel};
		$SEQ->{$gene}{lowq} = 0 if not defined $SEQ->{$gene}{lowq};
		$SEQ->{$gene}{pos} = 0 if not defined $SEQ->{$gene}{pos};
		$SEQ->{$gene}{neg} = 0 if not defined $SEQ->{$gene}{neg};
		my $gene2 = $SEQ->{$gene}{orig};
		my $text = "
- $gene (original name = $gene2):
Positive    = $SEQ->{$gene}{pos} (neg->pos = $SEQ->{$gene}{negpos})
Negative    = $SEQ->{$gene}{neg} (pos->neg = $SEQ->{$gene}{posneg})
UnkPos      = $SEQ->{$gene}{unkpos}
UnkNeg      = $SEQ->{$gene}{unkneg}
Used        = $SEQ->{$gene}{used}
Total       = $SEQ->{$gene}{total}
Too Short   = $SEQ->{$gene}{badlength}
High Indel  = $SEQ->{$gene}{indel}
Low Quality = $SEQ->{$gene}{lowq}
";
	#LOG($outLog, $text,"NA");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,header\tBAMFilename\tgene\ttotal\tused\tpositive\tnegative\tunkpos\tunkneg\ttooshort\tlowqual\n");
	#LOG($outReadLog, "footLoop.pl,read_passed_filter_gene,record\t$BAMFileName\t$gene\t$SEQ->{$gene}{total}\t$SEQ->{$gene}{used}\t$SEQ->{$gene}{pos}\t$SEQ->{$gene}{neg}\t$SEQ->{$gene}{unkpos}\t$SEQ->{$gene}{unkneg}\t$SEQ->{$gene}{badlength}\t$SEQ->{$gene}{lowq}\n\n");
	#LOG($outReadLog, "footLoop.pl,gene_skipped_low_read,$BAMFileName,$zero\n","NA");
	}

	$zero = $zero eq "" ? "(None)\n" : "\n$zero\n";

	#LOG($outLog, "\n(Genes that have <= 10 reads:) $LGN$zero$N\n","NA");
	LOG($outLog, date() . "\t${GN}SUCCESS$N: Total=$BAMStats->{total}, used=$BAMStats->{used}, Low Map Quality=$BAMStats->{lowq}, Too short=$BAMStats->{badlength}\n");
	#LOG($outLog, "Output: ${LGN}$filteredDir$N\n\n");

}

sub filter_BAM_read {
	my ($seq, $count, $readname, $gene, $reason) = @_;
	$count->{$reason} ++;
	$seq->{$gene}{$reason} ++;
	return($seq, $count);
}
