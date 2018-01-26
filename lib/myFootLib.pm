package myFootLib;

use strict; use warnings;
use vars qw(@EXPORT);
use parent "Exporter";

our $N		="\e[0m";
our $B		="\e[1m";
our $BL		="\e[0;30m";
our $BU		="\e[0;34m";
our $GN		="\e[0;32m";
our $YW		="\e[1;33m";
our $CY		="\e[0;36m";
our $RD		="\e[0;31m";
our $PR		="\e[0;35m";
our $BR		="\e[0;33m";
our $GR		="\e[0;37m";
our $WT		="\e[0m";
our $DGR		="\e[1;30m";
our $LBU		="\e[1;34m";
our $LGN		="\e[1;32m";
our $LCY		="\e[1;36m";
our $LRD		="\e[1;31m";
our $LPR		="\e[1;35m";
our $DIES   ="$LRD!!!\t$N";

our @EXPORT = qw(
colorseq
colorseqCG
color
colorCG
colorconv
revcomp
date
getuuid
getFullpathAll
myeval
checkBismarkIndex
getDate
getMD5
getFilename
getFullpath
parseExon
intersect
ex
def
makedir
linecount
mean
median
sd
se
tmm
tmmsd
tmmse
parse_cigar
LOG
makehash
DIE
$DIES
$N 
$B
$BL		 
$BU			 
$GN		 
$CY			 
$RD			 
$PR		 
$BR		 
$GR
$YW
$WT
$DGR
$LBU
$LGN
$LCY
$LRD
$LPR
);


#################################

sub DIE {
	return "\n$DIES Died at file $CY " . __FILE__ . " $N at line $LGN " . __LINE__ . " $N";
}

sub makehash {
   my ($keys, $vals) = @_;
   my $hash;
   for (my $i = 0; $i < @{$keys}; $i++) {
      $hash->{$keys->[$i]} = (defined $vals and defined $vals->[$i]) ? $vals->[$i] : 0;
   }
   return($hash);
}

sub getuuid {
	my ($uuid) = `uuidgen`;
	chomp($uuid);
	return $uuid;

}

sub checkBismarkIndex {
	my ($geneIndexes, $mainFolder, $outLog) = @_;
	my $bismark_folder = "$mainFolder/Bismark_indexes/footLoop/";
	mkdir "$mainFolder/" if not -d "$mainFolder/";
	mkdir "$mainFolder/Bismark_indexes" if not -d "$mainFolder/Bismark_indexes/";
	mkdir "$mainFolder/Bismark_indexes/footLoop/" if not -d "$mainFolder/Bismark_indexes/footLoop";
	mkdir $bismark_folder if not -d $bismark_folder;
	print $outLog "\n$LRD!!!$N FATAL ERROR: $geneIndexes is empty!\n" and die "\n" if -s $geneIndexes == 0;

	my ($md1) = `md5sum $geneIndexes` =~ /^(\w+)[ ]+/;
	$bismark_folder = "$mainFolder/Bismark_indexes/footLoop/$md1";
	mkdir $bismark_folder if not -d $bismark_folder;
	if (not -d $bismark_folder) {mkdir $bismark_folder; return(0, $bismark_folder);}
	if (-d $bismark_folder and not -e "$bismark_folder/Bisulfite_Genome/md5sum.txt") {mkdir $bismark_folder; return(0, $bismark_folder);}
	
	my @line = `cat $geneIndexes`;
	my $check = 0;
	for (my $i = 0; $i < @line; $i++) {
		chomp($line[$i]);
		my @arr = split("\t", $line[$i]);
		next if $line[$i] =~ /^track/;
		if (@arr >= 4 and $arr[2] =~ /^\d+$/ or $arr[3] =~ /^\d+/) {
			$check = 1;
		}
	}
	print $outLog "\n$LRD!!!$N FATAL ERROR: $geneIndexes doesn't seem to be a bed file (MUST contain$YW 4$N column of$YW tab$N separated:$LGN chromosome, start, end, and gene name$N!\n\n" and die "\n" if ($check == 0);
	return (1, $bismark_folder) if (-d $bismark_folder);
}
sub LOG {
   my ($outLog, $text, $STEP, $STEPCOUNT) = @_;
	if (defined $STEP) {
		$STEP += $STEPCOUNT if defined $STEPCOUNT;
		$text =~ s/\$STEP/$STEP/g;
	}
   print $outLog $text;
	print $text;
	return $STEP if defined $STEP;
}

sub getDate {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   my $date = "$mday $months[$mon] $year $hour:$min:$sec";
   my $timenow = $hour * 3600 + $min * 60 + $sec;
   return($date);
}


sub getFilename {
	my ($fh, $type) = @_;
	# Split folder and fullname, e.g. $fh = "folder1/folder2/filename.txt.bed"
	# @names    = (folder1,folder2,filename.txt.bed)
	# $fullname = "filename.txt.bed"
	# $folder   = "folder1/folder2"
	# @names2   = (filename, txt, bed)
	# shortname = "filename"
	$fh = "./" . $fh if $fh !~ /\//;
	if ($fh =~ /\/$/) {$fh =~ s/\/+$//;}
	my (@names)   = split("\/", $fh);
	my $fullname  = pop(@names);
	my $folder    = @names != 1 ? join("\/", @names) : $names[0] =~ /^\.\.$/ ? "../" : $names[0] eq "." ? "./" : $names[0] =~ /^\.\/$/ ? "./" : $names[0] =~ /^\w+/ ? "$names[0]/" : die "Can't extract folder from $fh :( (names = @names)\n";
	my @names2    = split(".", $fullname);
	my $shortname = @names2 == 0 ? $fullname : $names2[0];
#	print "Names=" . join(",", @names) . "\n\nFOLDER=$CY$folder$N, $shortname\n\n";
	return($shortname)                      if not defined($type);
	return($folder, $fullname)              if $type eq "folderfull";
	return($folder, $shortname)             if $type eq "folder";
	return($fullname)                       if $type eq "full";
	return($folder, $fullname, $shortname)  if $type eq "all";
}

sub getFullpath {
	my ($fh) = @_;
	
	my ($folder, $fullname) = getFilename($fh, "folderfull");
	my $currdir = `pwd`; chomp($currdir);
	die "Folder of fh=$fh (folder=$folder=) does not exist!\n" if not defined($folder) or not -d $folder;
	($folder) = `cd \"$folder\" && pwd`;
	chdir $currdir;
	chomp($folder);
	return("$folder/$fullname");
}

sub parseExon {
	my ($exonFile, $coor, $genewant, $outFolder, $length_seq) = @_;
	my ($genechr, $genebeg, $geneend) = split("\t", $coor); 
	my %exon;
	my %data;
	$genewant = uc($genewant);
	open (my $in, "<", $exonFile) or die;
	while (my $line = <$in>) {	
		chomp($line);
		my ($chr, $beg, $end, $gene, $num, $strand) = split("\t", $line);
		($chr, $end, $beg, $gene, $num, $strand) = split("\t", $line) if $beg > $end;
		$gene = uc($gene);
		$chr = "chr$chr" if $chr !~ /^chr/;
		my ($name, $number) = $gene =~ /^(.+)_(\d+)$/;
		$name = $gene if not defined($name);
		#print "$line\n" if $gene =~ /MYLIP/i;
		if ($chr eq $genechr and intersect($beg, $end, $genebeg, $geneend) == 1) {
			#print "$line\n";
			#next if uc($name) ne uc($genewant);
			die "Multiple exons at same start in gene $gene:\n$line\n\n" if defined($data{$name}{$num}{$beg});
			$data{$name}{$num}{$beg} = "$chr\t$beg\t$end\t$strand";
		}
	}
	close $in;
	open (my $in2, "<", $exonFile) or die;
	while (my $line = <$in2>) {	
		chomp($line);
		my ($chr, $beg, $end, $gene, $num, $strand) = split("\t", $line);
		$gene = uc($gene);
		$chr = "chr$chr" if $chr !~ /^chr/;
		my ($name, $number) = $gene =~ /^(.+)_(\d+)$/;
		$name = $gene if not defined($name);
		if (defined($data{$name}{$num})) {
			$data{$name}{$num}{$beg} = "$chr\t$beg\t$end\t$strand";
		}
	}
	close $in2;
	foreach my $name (keys %data) {
		foreach my $num (sort {$a <=> $b} keys %{$data{$name}}) {
			my %temp; my $count = 0;
			my $lastend = $genebeg;
			foreach my $beg (sort {$a <=> $b} keys %{$data{$name}{$num}}) {
				$count ++; 
				my ($chr, $beg2, $end, $strand) = split("\t", $data{$name}{$num}{$beg});
				#print "$name.$num:\t$chr\t$beg\t$end\t$strand\n";
				# Get intron if not exon 1
				if ($count > 1) {
					for (my $i = $lastend; $i < $beg; $i++) {
						$temp{$i} = $strand eq "+" ? 1 : -1;
					}
				}
				for (my $i = $beg; $i < $end; $i++) {
					$temp{$i} = $strand eq "+" ? 2 : -2;
				}
				$lastend = $end;
			}
			for (my $i = $genebeg; $i <= $geneend; $i++) {				
				my $value = defined($temp{$i}) ? $temp{$i} : 0; #center value
				my $value2 = $value == 2 ? 2 : $value == -2 ? -2 : 0; #up/down blocks value
				my $pos = $i - $genebeg;
				#print "$name $num POS $pos\tGENEBEG $genebeg\t$value\n" if $i - $genebeg < 10;
				$exon{$name}{$num}[0][$pos] = $value2;
				$exon{$name}{$num}[1][$pos] = $value;
			}
		}
	}
	my @final;
	my $genecount = 0;
	mkdir "$outFolder/exon/" if not -d "$outFolder/exon/";
	open (my $exonOut, ">", "$outFolder/exon/$genewant.exon");
	if (defined($exon{$genewant})) {
		foreach my $num (sort {$a <=> $b} keys %{$exon{$genewant}}) {
			print $exonOut "$genewant\t$genecount\t" . join("\t", @{$exon{$genewant}{$num}[0]}) . "\n"; $genecount ++;
			print $exonOut "$genewant\t$genecount\t" . join("\t", @{$exon{$genewant}{$num}[1]}) . "\n"; $genecount ++;
			print $exonOut "$genewant\t$genecount\t" . join("\t", @{$exon{$genewant}{$num}[0]}) . "\n"; $genecount ++;
			print $exonOut ".\t$genecount\t" . join("\t", (0) x scalar(@{$exon{$genewant}{$num}[0]})) . "\n"; $genecount ++;
		}
	}
	foreach my $name (keys %exon) {
		foreach my $num (sort {$a <=> $b} keys %{$exon{$name}}) {
			next if $name eq $genewant;
			print $exonOut "$name\t$genecount\t" . join("\t", @{$exon{$name}{$num}[0]}) . "\n"; $genecount ++;
			print $exonOut "$name\t$genecount\t" . join("\t", @{$exon{$name}{$num}[1]}) . "\n"; $genecount ++;
			print $exonOut "$name\t$genecount\t" . join("\t", @{$exon{$name}{$num}[0]}) . "\n"; $genecount ++;
			print $exonOut ".\t$genecount\t" . join("\t", (0) x scalar(@{$exon{$name}{$num}[0]})) . "\n"; $genecount ++;
		}
	}
	print $exonOut ".\t$genecount\t" . join("\t", (0) x $length_seq) . "\n"; $genecount ++;
	close $exonOut;
}

sub intersect {
	my ($start1, $end1, $start2, $end2) = @_;
	die "Died at intersect: start1 ($start1) can't be bigger than end1 ($end1)\n" if $start1 > $end1;
	die "Died at intersect: start2 ($start2) can't be bigger than end2 ($end2)\n" if $start2 > $end2;
	return(1) if ($start1 >= $start2 and $start1 <= $end2);
	return(1) if ($start2 >= $start1 and $start2 <= $end1);
	return(0);
}

sub linecount {
   my ($input1) = @_;
   my ($linecount) = `wc -l $input1` =~ /^(\d+)/ if $input1 !~ /.(bam|rmdup|gz|zip)$/;
      ($linecount) = `zcat $input1| wc -l` =~ /^(\d+)/ if $input1 =~ /.gz$/;
      ($linecount) = `samtools view $input1| wc -l` =~ /^(\d+)/ if $input1 =~ /.(bam|rmdup)$/;
   return $linecount;
}
sub makedir {
   my ($folder, $isFile) = @_;
	my $log = "";
   #$folder = getFullpath($folder);
   my @folders = split("/", $folder);
   die if @folders == 0;
   my $curr_folder = "";
   $log .= "FOLDER $LCY$folder$N is already exist!\n" and return $log if -e $folder;
   for (my $i = 0; $i < @folders; $i++) {
		last if $i == @folders-1 and defined $isFile;
      $curr_folder .= "$folders[$i]/";
      next if $curr_folder =~ /^\/$/;
      next if $curr_folder =~ /^(\/|\/home\/)$/;
      $log .= "$i: Undefined folder to add. Current folder=$LGN$curr_folder$N\n" and return $log if not defined $folders[$i];
      next if -d $curr_folder or -e $curr_folder;
      system("mkdir $curr_folder") == 0 or die "Failed to mkdir $curr_folder: $!\n";
      $log .= "$i. $curr_folder\n";
   }
   $log .= "FOLDER $LGN$curr_folder$N is made!\n" if -e $curr_folder;
   return ($log);
}

sub ex { # 0 means it's not exists (bad)
   my $files = $_[0];
   return 0 if not defined $files;
   if ($files =~ /^ARRAY/) {
      my @files = @{$files} if defined $files;
      return 0 if def(\@files) == 0;
      my @ex;
      for (my $i = 0; $i < @files; $i++) {
#        print "\t$i: $files[$i]\n";
         push(@ex, "$i. $files[$i]") if not -e $files[$i] and not -d $files[$i];
      }
   #  print "Files not exists:\n-" . join("\n-", @ex) . "\n" if @ex != 0 and defined $verb;
      return 1 if @ex == 0;
      return 0;
   }
   return 1;
}
sub def {
   my $files = $_[0];
   return 0 if not defined $files;
   if ($files =~ /^ARRAY/) {
      my @files = @{$files} if defined $files;
      my @def;
      for (my $i = 0; $i < @files; $i++) {
         push(@def, $i) if not defined $files[$i];
      }
      return 1 if @def == 0;
      return 0;
   }
   return 1;
}

sub median {
   my (@value) = @_;
   return "NA" if not defined($value[0]) or @value == 0;
   return(quartile(\@value, 0.5));
   #     return($value[0]) if @value == 1;
   #     @value = sort {$a <=> $b} (@value);
   #     my $median = 0.5*($value[int(@value/2)]+$value[int(@value/2)+1]) if @value % 2 == 0;
   #$median = $value[int(@value/2)+1] if @value % 2 != 0;
   #     return($median);
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
sub tmmse {
        my (@value) = @_;
   my $sd = sd(@value);
   my $perc = int(0.05*@value+0.5);
   my $total = @value - 2*$perc - 1;
   $total = 1 if $total == 0;
   my $se = $sd/sqrt($total);
   return($se);
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
sub se {
        my (@value) = @_;
   my $sd = sd(@value);
   my $total = @value;
   $total = 1 if $total == 0;
   my $se = $sd/sqrt($total);
   return($se);
}
sub sum {
   my (@value) = @_;
   return(0) if @value == 0;
   if (@value == 1 and $value[0] eq 'ARRAY') {
      @value = @{$value[0]};
   }
   my $sum = 0;
   for (my $i =0 ; $i < @value; $i++) {
      $sum += $value[$i];
   }
   return($sum);
}

sub parse_cigar {
   my ($cigar, $type) = @_;
   my @num = split("[A-Z]+", $cigar);
   my @alp = split("[0-9]+", $cigar);
   shift(@alp);
   my $length = 0;
   for (my $i = 0; $i < @num; $i++) {
      die "Died alp $i isn't I/D/M/N ($alp[$i]) at cigar:\n$cigar\n" if $alp[$i] !~ /^[IDMN]+$/;
      $length += $num[$i] if $alp[$i] ne "I";
   }
	return $length if defined $type and $type =~ /len/i;
	return(\@num, \@alp) if defined $type and ($type =~ /num/i or $type =~ /alp/i or $type =~ /cig/);
   return(\@num, \@alp, $length);
}

sub myeval {
	my ($var) = @_;
	my $count = 0;
	my $print = "";
	if ($var =~ /^ARRAY/) {
		for (my $i = 0; $i < @{$var}; $i++) {
			$count ++;
			my $val = defined ($var->[$i]) ? $var->[$i] : "VALUE_UNDEF";
			$print .= "myeval ARRAY \%var failed at i=$LGN$'i'$N, value='$val'\n" if not defined $var->[$i];
			$count ++ if not defined $var->[$i];
		}
	}
	elsif ($var =~ /^HASH/) {
		foreach my $key (sort keys %{$var}) {
			my $key2 = defined ($key) ? $key : "KEY_UNDEF";
			my $val = defined ($var->{$key}) ? $var->{$key} : "VALUE_UNDEF";
			$print .= "myeval HASH \%var failed at key='$key2', value='$val'\n" if not defined $key or not defined $var->{$key};
			$count ++ if not defined $key or not defined $var->{$key};
		}
	}
	else {
		my $val = defined ($var) ? $var : "VALUE_UNDEF";
		$print .= "myeval SCALAR \$var failed at value='$val'\n" if not defined $var;
		$count ++ if not defined $var;
	}
	return ($count, $print);
}

sub getFullpathAll {
   my @arr = @_;
   for (my $i = 0; $i < @arr; $i++) {
		my $file = $arr[$i];
      ($file) = getFullpath($file);
		$arr[$i] = $file;
   }
   return(@arr);
}

sub getMD5 {
   my ($file) = @_;
   my $filetemp = getFullpath($file);
   $filetemp = $file if not defined $filetemp or not -e $filetemp;
   if ($file !~ /md5(sum)?$/) {
      #print "FILE is not md5 ($file)\n";
      my ($folder, $filename) = getFilename($filetemp, "folder");
      $file = -e "$folder/.$filename.md5sum" ? "$folder/.$filename.md5sum" :
              -e "$folder/.$filename.md5"    ? "$folder/.$filename.md5" :
              -e "$folder/$filename.md5sum"  ? "$folder/$filename.md5sum" :
              -e "$folder/$filename.md5"     ? "$folder/$filename.md5" : "";
		if (not -e $filetemp) {
			die "file=\n$CY$filetemp$N\nHERE\n";
			return "NA";
		}
      system("md5sum $filetemp > $folder/.$filename.md5") == 0 or print date() . "mitochy::getMD5 $YW$filetemp$N: failed to md5sum $filetemp!\n" and return;
      $file = "$folder/.$filename.md5";
   }
	if (not -e $filetemp) {
		die "file=\n$CY$filetemp$N\nHERE\n";
		return "NA";
	}
   my @line = `cat $file`;
   return "NA" if @line == 0;
   return "NA" if $line[0] !~ /^[a-z0-9]+/i;
   chomp($line[0]);
   my ($md5, $file2) = split("[\t ]+", $line[0]);
   #print "\nMD5 = $md5\n\n";
   #print "\t$YW Warning$N: file=$CY$file$N, but file inside differs! ($LGN$file2$N)\n" if $file !~ /$file2/ and $file2 !~ /$file/;
   return($md5, $file2, $file);
}

sub date {
   my ($add, $color) = @_;
   ($color = $add and undef $add) if defined $add and $add =~ /^(color=|color|col=|col|y|c)/i;
   $add = 0 if not defined $add;
   return ("[$YW" . getDate($add) . "$N]: ") if not defined $color;
   return ("[" . getDate($add) . "]: ") if defined $color;
}

sub revcomp {
   my ($sequence) = @_;
   $sequence = uc($sequence);
   $sequence =~ tr/ATGC/TACG/;
   $sequence = reverse($sequence);
   return ($sequence);
}

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
sub makeopt {
	my 
(
'a' => $opt_a,
'b' => $opt_b,
'c' => $opt_c,
'd' => $opt_d,
'e' => $opt_e,
'f' => $opt_f,
'g' => $opt_g,
'h' => $opt_h,
'i' => $opt_i,
'j' => $opt_j,
'k' => $opt_k,
'l' => $opt_l,
'm' => $opt_m,
'n' => $opt_n,
'o' => $opt_o,
'p' => $opt_p,
'q' => $opt_q,
'r' => $opt_r,
's' => $opt_s,
't' => $opt_t,
'u' => $opt_u,
'v' => $opt_v,
'w' => $opt_w,
'x' => $opt_x,
'y' => $opt_y,
'z' => $opt_z,
'A' => $opt_A,
'B' => $opt_B,
'C' => $opt_C,
'D' => $opt_D,
'E' => $opt_E,
'F' => $opt_F,
'G' => $opt_G,
'H' => $opt_H,
'I' => $opt_I,
'J' => $opt_J,
'K' => $opt_K,
'L' => $opt_L,
'M' => $opt_M,
'N' => $opt_N,
'O' => $opt_O,
'P' => $opt_P,
'Q' => $opt_Q,
'R' => $opt_R,
'S' => $opt_S,
'T' => $opt_T,
'U' => $opt_U,
'V' => $opt_V,
'W' => $opt_W,
'X' => $opt_X,
'Y' => $opt_Y,
'Z' => $opt_Z,
'1' => $opt_1,
'2' => $opt_2,
'3' => $opt_3,
'4' => $opt_4,
'5' => $opt_5,
'6' => $opt_6,
'7' => $opt_7,
'8' => $opt_8,
'9' => $opt_9,
'0' => $opt_0
)

	return \%opt;
}
my ($useVars) = `grep 'use vars' $0` =~ /qw\((.+)\)/;
my @opts = split(" ", $useVars); @opts = sort @opts;
my %opts;
foreach my $opt (sort @opts) {
   my ($key) = $opt =~ /_(.)$/;
   $opts{$key} = $opt;
   print "$key -> $opts{$key}\n";
}
print "opts:\n" . join("\n", @opts) . "\n";die;
