#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_b $opt_L $opt_g $opt_i $opt_o $opt_B $opt_q $opt_x $opt_y $opt_r $opt_O);
getopts("vb:L:g:i:o:B:x:y:r:O:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib;
use FAlite;

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
LOG($outLog, "\n\n");
my $SEQ = parse_fasta($seqFile, $outLog, $minReadL);
parse_BAMFile($BAMFile, $seqFile, $filteredDir, $outLog);

LOG($outLog, "outDir=$outDir\n");
close $outLog;

sub parse_fasta {
   my ($seqFile, $outLog, $minReadL) = @_;
   open(my $SEQIN, "<", $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
   my $fasta = new FAlite($SEQIN);
   my $linecount = 0;
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
      $SEQ->{$gene}{minReadL}  = $minReadL; 
      $SEQ->{$gene}{total}     = 0;
      $SEQ->{$gene}{badlength} = 0;
      $SEQ->{$gene}{lowq}      = 0;
      $SEQ->{$gene}{used}      = 0;
      $SEQ->{$gene}{pos}       = 0;
      $SEQ->{$gene}{neg}       = 0;
      $SEQ->{$gene}{orig}      = $def;
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
	open($BAM, "samtools view $BAMFile|") or (LOG($outLog, "$LRD!!!$N\tFATAL ERROR: Could not open $BAMFile: $!") and exit 1) if $BAMFile =~ /.bam$/;

	## Some stats
	my $linecount   = 0;
	my $BAMStats = makehash(['total','used','diffgene','lowq','badlength']);

	while(my $line = <$BAM>) {
		$linecount ++;
		chomp($line);

		#####################
		# 1. Parse BAM line #
		#####################
		my @arr = split("\t", $line);
		LOG($outLog, "\t$YW$BAMFile$N: Done $linecount\n") if $linecount eq 1 or $linecount % 5000 == 0;

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
	LOG($outLog, "\n" . date() . "${LPR}parsing BAM$N: Parsed $LGN$linecount$N lines from bam file\n");
	LOG($outLog, "\n\tb.Logging BAM file $CY$BAMFile$N\n");
	foreach my $genez (sort keys %{$SEQ}) {
		my $outTXTFilePos  = "$filteredDir/$genez\_Pos.filtered.gz";
		my $outTXTFileNeg  = "$filteredDir/$genez\_Neg.filtered.gz";
		open (my $tempoutpos, ">", $outTXTFilePos) or die "Failed to write to \$outTXTFilePos $LCY$outTXTFilePos$N: $!\n";
		close $tempoutpos;
		open (my $tempoutneg, ">", $outTXTFileNeg) or die "Failed to write to \$outTXTFileNeg $LCY$outTXTFileNeg$N: $!\n";
		close $tempoutneg;
	}
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
	my ($BAMFileName) = getFilename($BAMFile);

	my $skipped = 0; my ($passedFilterP, $passedFilterN) = (0,0);
	open (my $inBAMFix, "zcat $filteredDir/$BAMFileName.fixed.gz|") or LOG($outLog, "Failed to open $filteredDir/$BAMFileName.fixed.gz: $!\n") and exit 1;
	my $linecount = 0;
	while (my $line = <$inBAMFix>) {
		chomp($line);
		LOG($outLog, date() . "${LPR}log_BamFile$N: Done $LGN$linecount$N\n") if $linecount % 100 == 0;
		$linecount ++;
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
		if ($strand eq 16) {
			my $outTXTFileNeg  = "$filteredDir/$genez\_Neg.filtered.gz";
			open (my $outTXT, "| gzip -f >> $outTXTFileNeg");
			print $outTXT "$read\tFN\t$pos\n";
			close $outTXT;
			$SEQ->{$genez}{neg} ++;
			$passedFilterN ++;
			$SEQ->{$genez}{posneg} ++ if $oldStrand eq 0;
		}
		if ($strand eq 0) {
			my $outTXTFilePos  = "$filteredDir/$genez\_Pos.filtered.gz";
			open (my $outTXT, "| gzip -f >> $outTXTFilePos");
			print $outTXT "$read\tFP\t$pos\n";
			close $outTXT;
			$SEQ->{$genez}{pos} ++;
			$passedFilterP ++;
			$SEQ->{$genez}{negpos} ++ if $oldStrand eq 16;
		}
		$BAMStats->{used} ++;
		$SEQ->{$genez}{used} ++;
	}
	close $inBAMFix;

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
		my $outTXTFilePos  = "$filteredDir/$gene\_Pos.filtered.gz"; 
		my $outTXTFileNeg  = "$filteredDir/$gene\_Neg.filtered.gz"; 

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
	}

	$zero = $zero eq "" ? "(None)\n" : "\n$zero\n";

	LOG($outLog, date() . "\t${GN}SUCCESS$N: Total=$BAMStats->{total}, used=$BAMStats->{used}, Low Map Quality=$BAMStats->{lowq}, Too short=$BAMStats->{badlength}\n");
}

sub filter_BAM_read {
	my ($seq, $count, $readname, $gene, $reason) = @_;
	$count->{$reason} ++;
	$seq->{$gene}{$reason} ++;
	return($seq, $count);
}
