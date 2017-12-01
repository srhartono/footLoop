package foot;

use strict; use warnings; 
#use Exporter;
#our @ISA = qw(Exporter);

use vars qw(@EXPORT);
use parent "Exporter";

############### COLOR #############
our $N      ="\e[0m";
our $B      ="\e[1m";
our $BL     ="\e[0;30m";
our $BU     ="\e[0;34m";
our $GN     ="\e[0;32m";
our $YW     ="\e[1;33m";
our $CY     ="\e[0;36m";
our $RD     ="\e[0;31m";
our $PR     ="\e[0;35m";
our $BR     ="\e[0;33m";
our $GR     ="\e[0;37m";
our $WT     ="\e[0m";
our $DGR    ="\e[1;30m";
our $LBU    ="\e[1;34m";
our $LGN    ="\e[1;32m";
our $LCY    ="\e[1;36m";
our $LRD    ="\e[1;31m";
our $LPR    ="\e[1;35m";
our $DIES   ="${LRD}!!!\t$N";

our @EXPORT = qw(getFilename getFullpath parseExon checkBismarkIndex $DIES $N $B $BL $BU $GN $CY $RD $PR $BR $GR $YW $WT $DGR $LBU $LGN $LCY $LRD $LPR);
#################################

sub checkBismarkIndex {
	my ($geneIndexes, $outLog) = @_;
	mkdir "/home/maligm/data/" if not -d "/home/maligm/data/";
	print "\n\nNot exist directory /home/maligm/data/Bismark_indexes/ \n\n" if not -d "/home/maligm/data/Bismark_indexes/";
	print "\n\nDirectory /home/maligm/data/Bismark_indexes/ EXIST!\n\n" if -d "/home/maligm/data/Bismark_indexes/";
	mkdir "/home/maligm/data/Bismark_indexes/" if not -d "/home/maligm/data/Bismark_indexes/";
	my $bismark_folder = "/home/maligm/data/Bismark_indexes/footLoop/";
	mkdir $bismark_folder if not -d $bismark_folder;
	print $outLog "\n$LRD!!!$N FATAL ERROR: $geneIndexes is empty!\n" and die "\n" if -s $geneIndexes == 0;

	my ($md1) = `md5sum $geneIndexes` =~ /^(\w+)[ ]+/;
	$bismark_folder = "/home/maligm/data/Bismark_indexes/footLoop/$md1";
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

	if ($fh =~ /\/$/) {$fh =~ s/\/+$//;}
	my (@names)   = split("\/", $fh);
	print "NAMES = @names\n";
	my $fullname  = pop(@names);
	my $folder    = @names == 1 ? "./" : join("\/", @names); 
	my @names2    = split(".", $fullname);
	my $shortname = $names2[0];
	
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
	chdir $folder if defined($folder) and -e $folder;
	($folder) = `pwd`;
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

