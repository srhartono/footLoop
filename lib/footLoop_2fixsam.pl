#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_n $opt_S $opt_c $opt_C $opt_o $opt_v $opt_n);
getopts("s:i:g:f:S:cCo:vn:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib;
use FAlite;

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

my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N $CY-n [folder of -n footLop.pl]$N $LGN-o$N [output dir]

";

my ($footLoop_folder) = $opt_n;
my ($footLoop_2fixsam_outDir) = $opt_o;
$footLoop_folder = "footLoop_folder_unknown" if not defined $opt_n;
my $footLoop_folder_forLog = $footLoop_folder;
$footLoop_folder_forLog =~ s/\/+/_/g;
$footLoop_folder_forLog =~ s/^\/+/SLASH_/;
#my $tempLog = "./.$footLoop_folder_forLog\_TEMPOUTLOG.txt";
#system("touch $tempLog") == 0 or print "Failed to write to $tempLog!\n";
#open (my $tempLogOut, ">", $tempLog) or print "Failed to write to $tempLog: $!\n";
#DIELOG($tempLogOut, "\nTEMPLOG:\n$tempLog\n$footLoop_folder: footLoop_2fixsam.pl: $usage") unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1 and -e $tempLog;
#DIELOG($tempLogOut, "\nTEMPLOG:\n$tempLog\n$footLoop_folder: footLoop_2fixsam.pl: please define output (-o)\n") if not defined $opt_o and -e $tempLog;
(print "\nfootLoop_2fixsam.pl: $usage\n" and exit) unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_n) == 1;
(print "\nfootLoop_2fixsam.pl: please define output (-o)\n" and exit) if not defined $opt_o;
makedir($footLoop_2fixsam_outDir);


###########
# LOGFILE #
###########

# log file
my $footLoop_logFile = "$footLoop_folder/logFile.txt";
my $footLoop_2fixsam_logFile = "$footLoop_folder/footLoop_2fixsam_logFile.txt";
open (my $outLog, ">", $footLoop_2fixsam_logFile) or die "Failed to write to footLoop_2fixsam_logFile: $!\n";
LOG($outLog, date() . "Logfile = $LRD$footLoop_logFile$N\n"); 

# parse footLoop logfile
my ($samFile, $seqFile, $genez) = parse_footLoop_logFile($footLoop_logFile, $outLog); 

# check sam file
LOG($outLog, date() . "Checking sam File =$LCY$samFile$N=\n");
check_file($samFile, "sam", $outLog); 

# check seq file
LOG($outLog, date() . "Checking seq File =$LCY$seqFile$N=\n");
check_file($seqFile, "seq", $outLog); 

# parse seq file and get chr etc
my %refs = %{parse_seqFile($seqFile)}; 
foreach my $chr (sort keys %refs) {
	$genez->{$chr} = @{$refs{$chr}};
}

my %out;
my %data; my $cons; my %strand;
my $linecount = 0;
my ($samFolder, $samName) = getFilename($samFile, "folderfull");
my $debugFile = "$footLoop_2fixsam_outDir/debug.txt";
my $outSam = "$footLoop_2fixsam_outDir/$samName.fixed";
open (my $outsam, ">", "$outSam") or die "Cannot write to $outSam: $!\n";
open (my $outdebug, ">", "$debugFile") or die "Cannot write to $debugFile: $!\n";
my ($total_read) = `awk '\$2 == 0|| \$2 == 16 {print}' $samFile | wc -l` =~ /^\s*(\d+)$/;
$linecount = 0;
my $in1;
open ($in1, "samtools view $samFile|") or die "Cannot read from $samFile: $!\n" if $samFile =~ /.bam$/;
open ($in1, "<", $samFile) or die "Cannot read from $samFile: $!\n" if $samFile =~ /.sam$/;
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	$linecount ++;
	LOG($outLog, date() . "\t$0: Parsed $LGN$linecount$N / $LCY$total_read$N\n","NA") if $linecount % 50 == 0;
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

	# 3. DETERMINING TYPE BASED ON CONVERSION
	my $type;

	# 3a 3_NONE: super low C->T conversion and G->A conversion then it's NONE
	if ($CT1 <= 5 and $GA1 <= 5) {
		$type = "3_NONE";
	}

	# 3b. 6_BOTH: otherwise, if C->T and G->A are within +/- 10% then of each other then it's BOTH
	elsif ($CT1 == $GA1 or ($CT1 > 5 and $GA1 > 5 and ($GA1 >= $CT1 * 0.9 and $GA1 <= $CT1 * 1.1) and ($CT1 >= $GA1 * 0.9 and $CT1 <= $GA1 * 1.1))) {
		$type = "6_BOTH";
	}
	# 3c. 1_SNEG and 6_SPOS: otherwise if one is strongly less than the other then STRONG "S" POS or NEG (SPOS or SNEG)
	# -> arbitrary criterias (CT vs GA for POS, and vice versa for NEG)
	#    1. CT 15+ vs GA 5-
	#    2. CT 5+ and ratio CT:GA is at least 3:1
	#    3. CT 20+ and ratio CT:GA is at least 2:1
	elsif (($CT1 > 15 and $GA1 < 5) or ($GA1 >= 5 and $CT1 / $GA1 > 3) or ($GA1 >= 20 and $CT1 / $GA1 >= 2)) {
		$type = "5_SPOS";
	}
	elsif (($GA1 > 15 and $CT1 < 5) or ($CT1 >= 5 and $GA1 / $CT1 > 3) or ($CT1 >= 20 and $GA1 / $CT1 >= 2)) {
		$type = "1_SNEG";
	}
	# 3d. 2_WNEG and 4_WPOS: otherwise just WEAK "W" POS or NEG (WPOS or WNEG)
	elsif ($CT1 < $GA1) {
		$type = "2_WNEG";
	}
	elsif ($CT1 > $GA1) {
		$type = "4_WPOS";
	}
	# 3e. 99_UNK: otherwise default is UNKNOWN
	else {
		$type = "99_UNK";
	}

	print $outsam "$read\t$type\t$strand\t$newstrand\t$chr\t$CTPrint\t$CT0,$CC0,$GA0,$GG0,$CT1,$CC1,$GA1,$GG1\n";
	LOG($outLog, date() . "file=$LCY$samFile$N, linecount=$linecount, read=$read, die coz no info\n") if not defined $GG1;

	## 3f. Below is for debug printing
	#print $outdebug ">$read,$type,OldStrand=$strand,NewStrand=$newstrand,$chr,CT0=$CT0,CC0=$CC0,GA0=$GA0,GG0=$GG0,CT1=$CT1,CC1=$CC1,GA1=$GA1,GG1=$GG1\n";
	#print $outdebug "$refPrint\n";
	#print $outdebug "$seqPrint\n";
	#print $outdebug "$CTPrint\n";
	#print $outdebug "$read\t$chr\tstrand=$strand, new=$newstrand\n" . join("", @{$ref3}) . "\n";
	#print $outdebug "CC = $CC1 / $CC0\n";
	#print $outdebug "CT = $CT1 / $CT0\n";
	#print $outdebug "GG = $GG1 / $GG0\n";
	#print $outdebug "GA = $GA1 / $GA0\n";
	#print $outdebug "REF: $refPrint\n";
	#print $outdebug "SEQ: $seqPrint\n";
	#print $outdebug "CON: " . join("", @{$CTcons}) . "\n";
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
			$tmm    = int(1000*tmm(@{$CT})+0.5)/1000;
			$mean   = int(1000*mean(@{$CT})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$CT})+0.5)/1000;
			$meanse = int(1000*se(@{$CT})+0.5)/1000;
		}
		print $outdebug "CT=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse, ";
		my $tot = $strand{$strand}{tot}{$type}; 
		($mean, $meanse, $tmm, $tmmse) = (0,0,0,0);
		if (defined $tot) {
			$tmm    = int(1000*tmm(@{$tot})+0.5)/1000;
			$mean   = int(1000*mean(@{$tot})+0.5)/1000;
			$tmmse  = int(1000*tmmse(@{$tot})+0.5)/1000;
			$meanse = int(1000*se(@{$tot})+0.5)/1000;
		}
		print $outdebug "tot=tmm=$tmm +/- $tmmse;mean=$mean +/- $meanse\n";
	}
}
exit 0;

# light quick dirty check if sam or seq file are sane
# sam is sane if there are at least 10 rows (or less if less than 20 reads) with more than 10 columns
# seq is sane if header is followed by seq and seq is ACTGUN (case-insensitive) at least 20 reads (or less if less than 20 reads in file)
sub check_file {
	my ($file, $type, $outLog) = @_;
	DIELOG($outLog, "footLoop_2fixsam.pl: $type file $file does not exist!\n") if ex($file) == 0;
	DIELOG($outLog, "footLoop_2fixsam.pl: $type file $file is empty!\n")       if -s $file  == 0;

	my $filetype = `file -b --mime-type $file`; chomp($filetype);
	my $cmd = ($file =~ /\.(rmdup|bam)$/ or $filetype =~ /(gzip|binary)/) ? "samtools view $file|" : "$file";
	my ($linecount, $check) = (0,0);
	my $currseq; #seq only

	my $total_line = 0;
	my $checkfileIn;
	if ($file =~ /\.(rmdup|bam)$/) {
		($total_line) = `samtools view $file| wc -l` =~ /^\s*(\d+)/;
		LOG($outLog, "samtools view $file| wc -l = $total_line\n","NA");
		open ($checkfileIn, "samtools view $file|") or DIELOG($outLog, "footLoop_2fixsam.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n");
	}
	elsif ($file =~ /\.gz$/ or ($filetype =~ /(gzip|binary)/ and $file !~ /\.(rmdup|bam)$/)) {
		($total_line) = `zcat < $file| wc -l` =~ /^\s*(\d+)/;
		open ($checkfileIn, "zcat < $file|") or DIELOG($outLog, "footLoop_2fixsam.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n");
	}
	else {
		($total_line) = `wc -l $file` =~ /^\s*(\d+)/;
		open ($checkfileIn, "<", $file) or DIELOG($outLog, "footLoop_2fixsam.pl: Failed to read from filetype=$filetype, file=$LCY$file$N: $!\n") 
	}
	
	while (my $line = <$checkfileIn>) {
		$linecount ++;
		chomp($line); my @arr = split("\t", $line);
		last if $check >= 20 or $linecount >= $total_line;
		if ($type eq "sam") {
			if (@arr > 10) {
				$check = $check == 0 ? 2 : $check + 1;
			}
		}
		elsif ($type eq "seq") {
			LOG($outLog, "linecount = $linecount check=$check line = $line\n","NA");
			while ($line !~ /^>/) {
				last if $check >= 20 or $linecount >= $total_line;
				LOG($outLog, "\tlinecount=$linecount, check=$check\n", "NA");
				$currseq .= $line;
				$line = <$checkfileIn>; 
				chomp($line); $linecount ++;
			}
			last if $check >= 20 or $linecount >= $total_line;
			if (defined $currseq and $currseq =~ /^[ACTGUN]+$/i) {
				LOG($outLog, __LINE__ . "before: linecount=$linecount check=$check file=$file type=$type\n","NA");
				$check = $check == 0 ? 0 : $check + 1;
				LOG($outLog, __LINE__ . "after: linecount=$linecount check=$check\n","NA");
				DIELOG($outLog, __LINE__ .  "footLoop_2fixsam.pl: linecount=$linecount file=$file type=$type, check=$check, check_file failed (does not seem to be a seq file\n\t-> fasta file start with non-header!\n\n") if $check == 0;
				DIELOG($outLog, __LINE__ .  "footLoop_2fixsam.pl: linecount=$linecount file=$file type=$type check=$check line=$line corruped fasta file!\n\t-> multiple headers in a row\n\n") if $check % 2 != 0;
				undef $currseq;
			}
			if ($linecount < $total_line and defined $line and $line =~ /^>/) {
				LOG($outLog, __LINE__ . "linecount=$linecount check=$check file=$file type=$type\n","NA");
				$check ++;
				# line will always be parsed header (1,3,5,etc) then seq (2,4,6,etc) so if header is even number then corrupted fasta file
				DIELOG($outLog, __LINE__ .  "footLoop_2fixsam.pl: linecount=$linecount file=$file type=$type check=$check line=$line corruped fasta file!\n") if $check % 2 == 0;
			}
		}
		last if $check >= 20 or $linecount >= $total_line;
	}
	DIELOG($outLog, __LINE__ .  "footLoop_2fixsam.pl: linecount=$linecount file=$file type=$type, check=$check, check_file failed (does not seem to be a sam file\n\t-> there's no read with more than 10 collumns!\n") if $check == 0;
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
		#LOG($outLog, date() . "REF $def length=" . length($seq) . " SEQ = @seq\n\n") if $def eq "CALM3";
   }
   close $in;
	return(\%ref);
}

sub parse_footLoop_logFile {
	my ($footLoop_logFile) = @_;
	DIELOG($outLog, "footLoop_2fixsam.pl: \n\nCan't find $footLoop_logFile! Please run footLoop.pl first before running this!\n\n") if not -e $footLoop_logFile;
#ADD
	$footLoop_logFile = "$footLoop_folder/.PARAMS";
	LOG($outLog, date() . "\t\tLOGFILE=$footLoop_logFile\n");
	my ($samFile, $seqFile, $genez);
	if (-e $footLoop_logFile) {
		my @line = `cat $footLoop_logFile`;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i]; chomp($line);
			if ($line =~ /footLoop.pl,samFile/) {
				($samFile) = $line =~ /samFile,(.+)$/;
			}
			if ($line =~ /footLoop.pl,seqFile/) {
				($seqFile) = $line =~ /seqFile,(.+)$/;
			}
			if ($line =~ /footLoop.pl,geneIndexFile/) {
				my ($geneIndexFile) = $line =~ /geneIndexFile,(.+)$/;
				my @line = `cut -f1-4 $geneIndexFile`;
				foreach my $line (@line) {
					chomp($line); next if $line =~ /^#/;
					my ($chr, $beg, $end, $gene) = split("\t", $line);
					my $length = $end - $beg;
					$gene = uc($gene);
					$genez->{$gene} = $length;
					LOG($outLog, date() . "\t\tfootLoop_2fixsam.pl: gene $gene = $length bp\n");
				}
			}
			
#			if ($line =~ /^!samFile=/) {
#				($samFile) = $line =~ /^!samFile=(.+)$/;
#			}
#			if ($line =~ /^!seqFile=/) {
#				($seqFile) = $line =~ /^!seqFile=(.+)$/;
#			}
#			if ($line =~ /gene=.+length=/) {
#				my ($gene, $length) = $line =~ /^.+gene=(.+) length=(\d+)$/;
#				die if not defined $gene or not defined $length;
#				$genez->{$gene} = $length;
#				print "gene $gene = $length bp\n";
#			}
#			last if $line =~ /footLoop_2fixsam.pl/;
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
			DIELOG($outLog, "footLoop_2fixsam.pl: Undefined ref2 at j=$j\n") if not defined $ref2->[$j];
			$reftemp .= $ref2->[$j];
			DIELOG($outLog, "footLoop_2fixsam.pl: Undefined seq2 at j=$j\n") if not defined $seq2->[$j];
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
		($ref) = $alp[$i] eq "I" ? ins($ref, $refpos, "-", $num[$i], "ref") : $ref;
		($seq) = $alp[$i] eq "D" ? ins($seq, $seqpos, "-", $num[$i], "seq") : $seq; 
		$refpos += $num[$i];
		$seqpos += $num[$i];
		$insref += $num[$i] if $alp[$i] eq "I";
		$insseq += $num[$i] if $alp[$i] eq "D";
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
	DIELOG($outLog, "footLoop_2fixsam.pl: X=$X, H=$H, Z=$Z, U=$u, x=$x, h=$h, z=$z, u=$u, dot=$dot, length=$length, length2=$length2\n\n") if $length ne $length2;
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
		DIELOG($outLog, "footLoop_2fixsam.pl: Can't find gene $chr in $seqFile!\n") if not defined $refs{$chr};
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

