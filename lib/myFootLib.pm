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
parse_readName
parse_footPeak_logFile
parse_indexFile
myformat
perc
fold
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
getMD5_simple
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
parseName
LOG
makehash
DIE
DIELOG
makeOutDir
getFlag
shuffle
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

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";

#################################

sub get_pcb_readname {
	my ($outLog) = @_;
	my %data;
	LOG($outLog, "\n\n$0::get_pcb_readname\n\n") if defined $outLog;
	LOG($outLog, date() . "\ncurl https://s3-us-west-1.amazonaws.com/muhucsc/pcb_readname.tsv\n\n") if defined $outLog;
	my @line = `curl https://s3-us-west-1.amazonaws.com/muhucsc/pcb_readname.tsv 2>&1`;
	foreach my $line (@line) {
		next if $line !~ /^PCB/;
		my ($PCB, $readName) = split("\t", $line);
		my ($id1a, $id1b) = $readName =~ /^.*m(\d+)_(\d+)_\d+/;
		DIELOG($outLog, "Failed to parse id1/2/3 from readName=$LCY$readName$N\n") if (not defined $id1a or not defined $id1b) and defined $outLog;
		die "Failed to parse id from readName=$LCY$readName$N\n" if (not defined $id1a or not defined $id1b) and not defined $outLog;
		my $id1 = $id1a . $id1b;
		my $id = "$id1"; $id =~ s/_//g;
		LOG($outLog, date() . "$0::get_pcb_readname: $PCB\t$id\n","NA") if defined $outLog;
		$data{readName}{$id} = $PCB;
		$data{PCB}{$PCB} = $id;
	}
	LOG($outLog, date() . "$0::get_pcb_readname: Parsed $LGN" . scalar(keys %{$data{PCB}}) . "$N PCBs from pcb_readname.tsv\n") if defined $outLog;
	return \%data;
}
sub parse_readName_from_fastq {
	my ($fastqFile, $outLog) = @_;
	my %name;
	my $in;
	open ($in, "<", $fastqFile) or DIELOG($outLog, "Failed to read from $fastqFile: $!\n") if defined $outLog;
	open ($in, "<", $fastqFile) or die "Failed to read from $fastqFile: $!\n" if not defined $outLog;
	while (my $line = <$in>) {
		chomp($line);
		my $def = $line;
		$line = <$in>; chomp($line);
		my $seq = $line;
		$line = <$in>; chomp($line);
		my $plus = $line;
		$line = <$in>; chomp($line);
		my $qua = $line;
		my ($id, $id1, $id2, $id3) = parse_readName($def);
		$name{$id1} ++;
	}
	LOG($outLog, "#id\tfreq\n") if defined $outLog;
	print "#id\tfreq\n" if not defined $outLog;
	foreach my $id (sort {$name{$b} <=> $name{$a}} keys %name) {
		LOG($outLog, "$id\t$name{$id}\n") if defined $outLog;
		print "$id\t$name{$id}\n" if not defined $outLog;
	}
	return \%name;
}

sub parse_readName {
	my ($readName, $outLog) = @_;
	DIELOG($outLog, "readName not defined!\n") if defined $outLog and not defined $readName;
	die "readName not defined!\n" if not defined $outLog and not defined $readName;
	return -1 if not defined $readName;
	my ($id1a, $id1b, $junk, $id2, $id3) = $readName =~ /^.*m(\d+)_(\d+)_\d+(_c\w+_\w+_\w+)?\/(\d+)\/(ccs|\d+_\d+)/;
	   ($id1a, $id1b, $junk, $id2, $id3) = $readName =~ /^.*m(\d+)_(\d+)(_\w+)?\/(\d+)\/(ccs|\d+_\d+)/ if not defined $id1a;
	   ($id1a, $id1b, $junk, $id2, $id3) = $readName =~ /^.*m(\d+)_(\d+)(_\w+)?\/(\d+)\/(ccs|\d+)/ if not defined $id1a;
	$id3 = (not defined $id3) ? 0 : $id3 eq "ccs" ? 0 : $id3;
	DIELOG($outLog, "Failed to parse id1/2/3 from readName=$LCY$readName$N\n") if (not defined $id1a or not defined $id1b or not defined $id2 or not defined $id3) and defined $outLog;
	my $id1 = $id1a . $id1b;
	my $id = "$id1$id2$id3"; $id =~ s/_//g;
	return ($id, $id1, $id2, $id3);

}
sub shuffle {
        my ($value, $times) = @_;
   my @value = @{$value};
   $times = @value if not defined($times) or $times !~ /^\d+$/;
        #print "Before: @value\n";
        for (my $i = 0; $i < $times; $i++) {
                my $rand1 = int(rand(@value));
                my $rand2 = int(rand(@value));
                my $val1 = $value[$rand1];
                my $val2 = $value[$rand2];
                $value[$rand1] = $val2;
                $value[$rand2] = $val1;
        }
        #print "After: @value\n";
        return(@value);
}

sub prettyPrint {
	my ($texts) = @_;
	my %len;
	my @texts = split("\n", $texts); #row
	for (my $i = 0; $i < @texts; $i++) { 
		chomp($texts[$i]);
		my @text = split("\t", $texts[$i]);#col
		for (my $j = 0; $j < @text; $j++) { 
			my $len = length($text[$j]);
			$len{$j} = $len if not defined $len{$j};
			$len{$j} = $len if $len{$j} < $len;
		}
	}

	my $newtexts;
	for (my $i = 0; $i < @texts; $i++) { 
		my @text = split("\t", $texts[$i]);#col
		for (my $j = 0; $j < @text; $j++) {
			my $len = $len{$j}; 
#			print "j = $j len = $len text=$text[$j]\n";
			my $newtext = $text[$j] . join("", ((" ") x ($len-length($text[$j])) ));
			$newtexts .= $j == @text - 1 ? "$newtext\n" : "$newtext\t";
		}
	}
	return($newtexts);
}

sub defFlag {
	my @PEAK = ("PEAK","NOPK");
	my @TEMP = ("", "_TEMP");
	my @RCONV = ("", "_RCONV");
	my @CPG  = ("", "_C");
	my @ALL = ("ALL");
	return(\@PEAK, \@TEMP, \@RCONV, \@CPG, \@ALL);
}
sub makeOutDir {
	my ($resDir, $bool) = @_;
	my ($PEAK, $TEMP, $RCONV, $CPG, $ALL) = defFlag();
	$bool = 0 if not defined $bool;
	$resDir =~ s/\/+$//;
	my $OUTDIRS;
	foreach my $PEAKS ( @{$PEAK}[0..(@{$PEAK}-1)] ) {
		foreach my $TEMPS ( @{$TEMP}[0..(@{$TEMP}-1)] ) {
			foreach my $RCONVS (@{$RCONV}[0..(@{$RCONV}-1)]) {
				foreach my $CPGS (@{$CPG}[0..(@{$CPG}-1)]) {
					my $outDir = $resDir . "/$PEAKS$TEMPS$RCONVS$CPGS";
					my $outDirName = "/$PEAKS$TEMPS$RCONVS$CPGS";
					$OUTDIRS->{$outDirName} = $outDir;
					makedir($outDir) if not -d $outDir and $bool eq 0;
					makedir("$outDir/ALL") if not -d "$outDir/ALL" and $bool eq 0;
				}
			}
		}
	}
	$OUTDIRS->{ALL} = "$resDir/$ALL->[0]";
	makedir("$resDir/$ALL->[0]/") if not -d "$resDir/$ALL->[0]/" and $bool eq 0;
	makedir("$resDir/$ALL->[0]/ALL/") if not -d "$resDir/$ALL->[0]/ALL/" and $bool eq 0;
	return($OUTDIRS);
}

sub getFlag {
	my ($file, $geneStrand, $readStrand, $rconvType) = @_;
	my ($PEAK, $TEMP, $RCONV, $CPG, $ALL) = defFlag();
	my $flag = $readStrand eq "Unk" ? "UNK" : $file =~ /NOPK/ ? "NOPK" : "PEAK";
	return ($ALL->[0]) if $readStrand eq "Unk";
	my $temp  = $geneStrand eq $readStrand ? $TEMP->[0] : $TEMP->[1];	
	my $rconv = (($rconvType =~ /^C/ and $readStrand eq "Pos") or ($rconvType =~ /^G/ and $readStrand eq "Neg")) ? $RCONV->[0] : $RCONV->[1];
	my $cpg   = $rconvType !~ /^(CG|GC)$/ ? $CPG->[0] : $CPG->[1];
	return($flag . $temp . $rconv . $cpg);
}
sub parse_indexFile {
	my ($indexFile) = @_;
	my %data;
	open (my $in, "<", $indexFile) or die "Failed to read from $indexFile: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /^\#/;
		my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
		die "\n\nmyFootLib.pm parse_indexFile: ERROR: chr undefined from indexFile $LCY$indexFile$N line=\n\n$line\n\n" if not defined $chr;
		die "\n\nmyFootLib.pm parse_indexFile: ERROR: beg undefined from indexFile $LCY$indexFile$N line=\n\n$line\n\n" if not defined $beg;
		die "\n\nmyFootLib.pm parse_indexFile: ERROR: end undefined from indexFile $LCY$indexFile$N line=\n\n$line\n\n" if not defined $end;
		die "\n\nmyFootLib.pm parse_indexFile: ERROR: gene undefined from indexFile $LCY$indexFile$N line=\n\n$line\n\n" if not defined $gene;
		die "\n\nmyFootLib.pm parse_indexFile: ERROR: strand undefined from indexFile $LCY$indexFile$N line=\n\n$line\n\n" if not defined $strand;
		%{$data{$gene}} = ("chr"=>$chr,"beg"=>$beg,"end"=>$end,"strand"=>$strand);
	}
	close $in;
	return \%data;
}

sub parseName {
	my ($filename, @info) = @_;
	$filename =~ s/\/+/\//g;
	if ($filename =~ /\//) {
		my @filename = split("/", $filename);
		$filename = pop(@filename);
	}
	my ($label, $gene, $strand, $window, $thres, $type) = $filename =~ /^(.+)_gene(.+)_(Pos|Neg|Unk)_(.+)_(.+)_(CH|CG|GH|GC)/;
	my ($pcb, $bc, $plasmid, $desc) = ("", "", "", "");
	if ($label =~ /^(.+)_bc.+_plasmid.+_desc.+$/) {
		($pcb, $bc, $plasmid, $desc) = $label =~ /^(.+)_bc(.+)_plasmid(.+)_desc(.+)$/;
		die "\n\nmyFootLib::parseName: filename=$LGN$filename$N.\n\nCannot parse bc, plasmid, desc from label=$LPR$label$N\n\n" if not defined $bc or not defined $plasmid or not defined $desc;
	}
	if (not defined $label or not defined $gene or not defined $strand or not defined $window or not defined $thres or not defined $type) {
		print "Cannot parse label gene strand window thres type from filename=$LCY$filename$N\n\nMake sure that filename format is: (.+)_gene(.+)_(Pos|Neg|Unk)_(.+)_(.+)_(CH|CG|GH|GC)\n\n";
		return -1;
	}
	$bc = "" if not defined $bc;
	$plasmid = "" if not defined $plasmid;
	$desc = "" if not defined $desc;
	$pcb = "" if not defined $pcb;
	$pcb = $label if not defined $pcb or $pcb eq "";
	my @data = ($label, $gene, $strand, $window, $thres, $type, $bc, $plasmid, $desc, $pcb);
	my %data = (
		"label" => $label,
		"gene" => $gene,
		"strand" => $strand,
		"window" => $window,
		"thres" => $thres,
		"type" => $type,
		"pcb" => $pcb,
		"bc" => $bc,
		"desc" => $desc,
		"plasmid" => $plasmid,
		"array" => \@data
	);
	return (\%data);
}
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

	my $tomd5 = $geneIndexes;
	my $md5res = "";
   if ($md5script eq "md5sum") {
      $md5res = `$md5script $tomd5`;
		die "Failed to $md5script $tomd5: $!\n" if not defined $md5res;
   }
   if ($md5script eq "md5") {
      ($md5res) = `$md5script $tomd5` =~ /^.+\= (\w+)$/; 
		die "Failed to $md5script $tomd5: $!\n" if not defined $md5res;
   }
	my $md1 = $md5res;
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
	if (defined $STEP and $STEP eq "NA") {
		print $outLog $text;
		return;
	}
	if (defined $STEP) {
		$STEP += $STEPCOUNT if defined $STEPCOUNT;
		$text =~ s/\$STEP/$STEP/g;
	}
	if (defined $outLog) {
		$text = "NO TEXT!?" if not defined $text;
	   print $outLog $text;
	}
	else {
		print "$LRD\tmyFootLib.pm:: LOG \$outLog isn't defined!! Only printing text$N\n";
	}
	print $text;
	return $STEP if defined $STEP;
}

sub DIELOG {
   my ($outLog, $text) = @_;
	#my ($scriptFileName) = getFilename($script, "fullname");
#	foreach my $args (@args) {#
		
#	}
	my ($STEP, $STEPCOUNT);
	if (defined $STEP) {
		$STEP += $STEPCOUNT if defined $STEPCOUNT;
		$text =~ s/\$STEP/$STEP/g;
	}
	print $outLog $text;
	print $text;
	die "DEAD\n";
}

sub getDate {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   my $date = "$mday $months[$mon] $year $hour:$min:$sec";
   my $timenow = $hour * 3600 + $min * 60 + $sec;
   return($date);
}

sub getDate_simple {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	($year) = $year =~ /^\d\d(\d+)$/;
	($mon) = $mon + 1 < 10 ? 0 . ($mon+1) : $mon + 1;
	return $year . $mon . $mday;
 #  my $date = "$mday $months[$mon] $year $hour:$min:$sec";
  # my $timenow = $hour * 3600 + $min * 60 + $sec;
  # return($date);
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
	$fh = "./" . $fh if $fh !~ /\//;
	my (@names)   = split("\/", $fh);
	my $fullname  = pop(@names);
	my $folder    = @names != 1 ? join("\/", @names) : $names[0] =~ /^\.\.$/ ? "../" : $names[0] eq "." ? "./" : $names[0] =~ /^\.\/$/ ? "./" : $names[0] =~ /^\w+/ ? "$names[0]/" : die "Can't extract folder from $fh :( (names = @names)\n";
	my @names2    = split(".", $fullname);
	my $shortname = @names2 == 0 ? $fullname : $names2[0];
#	print "Names=" . join(",", @names) . "\n\nFOLDER=$CY$folder$N, $shortname\n\n";
	$folder = "./" if $folder eq "";
	return($shortname)                      if not defined($type);
	return($folder, $fullname)              if $type eq "folderfull";
	return($folder, $shortname)             if $type eq "folder";
	return($fullname)                       if $type eq "full" or $type eq "fullname";
	return($folder, $fullname, $shortname)  if $type eq "all";
	print "$0::getFilename: fh=$fh type=$type TYPE DOES NOT EXIST!\n";
	return $shortname;
}

sub getFullpath {
	my ($fh, $die) = @_;
	$fh =~ s/\/+/\//g;
	$fh = "./" if $fh eq "." or $fh eq "." or $fh eq "" or not defined $fh;
	if (-d $fh) {
		my $folder = `cd $fh && pwd`;
		chomp($folder);
		die "fh = $fh, folder = $folder\n" if defined $die;
		return $folder;
	}
	$fh = "./$fh" if $fh !~ /^(\/|\.)/;
	my ($folder, $fullname) = getFilename($fh, "folderfull");
	my $currdir = `pwd`; chomp($currdir);
	die "Folder of fh=$fh (folder=$folder=) does not exist!\n" if not defined($folder) or not -d $folder;
	($folder) = `cd \"$folder\" && pwd`;
	chdir $currdir;
	chomp($folder);
	$folder = "./" if $folder eq "";
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
   my ($linecount) = `wc -l $input1` =~ /^\s*(\d+)/ if $input1 !~ /.(bam|rmdup|gz|zip)$/;
      ($linecount) = `zcat < $input1| wc -l` =~ /^\s*(\d+)/ if $input1 =~ /.gz$/;
      ($linecount) = `samtools view $input1| wc -l` =~ /^\s*(\d+)/ if $input1 =~ /.(bam|rmdup)$/;
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

sub quartile {
   my ($value, $p, $sort) = @_;
   die "${YW}mitochy::quartile$N: Probability Quartile has to be between 0 to 1! You submitted:$LRD $p$N\n" unless ($p >= 0 and $p <= 1);
   my @value = @{$_[0]}; my $total = @value;
#  print "prob = $p, total = " . scalar(@value) . "\n";
   @value = sort {$a <=> $b} (@value) if not defined $sort;
   my $i = (@value+1) * ($p) - 1;
   my $j = int($i);
#  print "\ti = $i = {($total+1) * $p - 1)}\n\tj = $j = int of i=$i\n";
#  print "\tvalue #j=$j * (i=$i - j=$j) $LGN PLUS $N value int($i+1) * ( 1 - (i=$i - j=$j))\n";
#  print join("\n",@value) . "\n";
   my $temp = int($i+1) > @value-1 ? @value-1 : int($i+1); $temp = 0 if $temp < 0;
   die "array of values does not have any value!\n" if @value == 0;
   die "total value is " . scalar(@value) . " but value $temp isn't defined!\n" if not defined $value[$temp];
## print "Val1 = value[j]=$value[$j] * (i=$i - j=$j)\n";
#  print "Val2 = value[temp]=$value[$temp] * (1 - (i=$i - j=$j))\n";
   my $val1 = $value[$j] * ($i-$j);
#  my $val2 = $value[int($i+1)] * (1 - ($i - $j));
   my $val2 = $value[$temp] * (1 - ($i - $j));
#  print "Val1 = value[j]=$value[$j] * (i=$i - j=$j) = $val1\n";
#  print "Val2 = value[temp]=$value[$temp] * (1 - ($i - j=$j)) = $val2\n";
#  print "\n\n";
   my $j2 = $j > @value-1 ? @value-1 : $j;
   my $i2 = $i > @value-1 ? @value-1 : $i;
   return($value[$j2] * ($i2 - $j2) + $value[$temp] * (1 - ($i2 - $j2)) );
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

sub max {
   my (@value) = @_;
   return(0) if @value == 0;
   if (@value == 1 and $value[0] eq 'ARRAY') {
      @value = sort {$b <=> $a} @value;
		return $value[0];
   }
}
sub min {
   my (@value) = @_;
   return(0) if @value == 0;
   if (@value == 1 and $value[0] eq 'ARRAY') {
      @value = sort {$a <=> $b} @value;
		return $value[0];
   }
}

sub myformat {
	my ($number, $note) = @_;
	$note = "" if not defined $note;
	$number =~ s/^\\\-/-/;
	print "\n\n$LRD WARNING!$N at myFOotLib.pl::myformat number=$CY$number$N does not look like number!\nNOTE: \"$note\"\n\n" and return ($number) if $number !~ /^\-?\d+\.?\d*e?\-?\d*$/;
	return $number if $number == 0;
	my $mod = $number < 0 ? -1 : 1;
	$number *= -1 if $number < 0;
	if ($number > 5) {
		return($mod*int($number*10+0.5)/10);
	}
	elsif ($number > 0.1 and $number <= 5) {
		return($mod*int($number*100+0.5)/100);
	}
	elsif ($number > 0.0000001 and $number !~ /e/i) {
		return($mod*abs($number));
#		my ($zero) = $number =~ /0\.(0+)[1-9]/;
#		($number) = $number =~ /^(0\.0+[1-9]\d?)/;
#		$zero = 10**(length($zero) + 3);
#		return($mod*int($number*10**$zero +0.5)/10**$zero);
	}
	else {# ($number =~ /^\-?\d+\.?\d*e\-?\d+/) 
		return ($mod*abs($number));
	}
}
sub perc {
	my ($peak, $total) = @_;
	return(int($peak/$total*1000+0.5)/10);
}
sub fold {
	my ($shuf, $orig) = @_;
	return(int(($shuf+5)/($orig+5)*1000+0.5)/1000);
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
			$print .= "myeval ARRAY \%var failed at i=$LGN'$i'$N, value='$val'\n" if not defined $var->[$i];
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

sub getMD5_simple {
	my ($file) = @_;
	my ($md5) = `$md5script $file` =~ /^(\w+)\s+/;
	return ($md5);
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
			print "$CY$filetemp$N\nfile does not exist!\n\n";
			return "NA";
		}
		my $tomd5 = $filetemp;
		my $md5res = "";
	   if ($md5script eq "md5sum") {
	      $md5res = `$md5script $tomd5`;
			die "Failed to $md5script $tomd5: $!\n" if not defined $md5res;
	   }
	   if ($md5script eq "md5") {
	      ($md5res) = `$md5script $tomd5` =~ /^.+\= (\w+)$/; 
			die "Failed to $md5script $tomd5: $!\n" if not defined $md5res;
	   }
      system("echo '$md5res\t$tomd5' > $folder/.$filename.md5") == 0 or print date() . "mitochy::getMD5 $YW$filetemp$N: failed to $md5script $filetemp!\n" and return;

      $file = "$folder/.$filename.md5";
   }
	if (not -e $filetemp) {
		print "$CY$filetemp$N\nfile does not exist!\n\n";
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

sub parse_footPeak_logFile {
	my ($footPeak_logFile, $footPeak_Folder, $genewant, $outLog) = @_;
   my %coor;
   my @lines   = `cat $footPeak_logFile`;
   LOG($outLog, "$YW$0$N ::parse_footPeak_logFile: Parsing footpeak logfile $footPeak_logFile\n","NA");
   my ($label) = `cat $footPeak_Folder/.LABEL`; chomp($label);
   DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_logFile!\n\n") if not -e $footPeak_logFile;
   DIELOG($outLog, "\n\ndied at footPeak_graph.pl: can't find $footPeak_Folder/.LABEL!\n\n") if not -e "$footPeak_Folder/.LABEL";
   my ($thres, $window);
   foreach my $line (@lines) {
      chomp($line);
      if ($line =~ /^[ \t]*def=.+, coor=.+/) {
         $line =~ s/(^\s+|\s+$)//g;
         my ($gene, $CHR, $BEG, $END, $GENE, $VAL, $STRAND) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\-?\d+\.?\d*), ([\+\-])$/;
         LOG($outLog, "gene=$LCY$gene$N,chr=$LCY$CHR$N,beg=$LCY$BEG$N,end=$LCY$END$N,gene=$LCY$GENE$N,val=$LCY$VAL$N,strand=$LCY$STRAND\n","NA");
         if (defined $genewant and $gene !~ /$genewant/i) {
            LOG($outLog, date() . " Skipped $LCY$gene$N as it doesn't contain $LGN-G $genewant$N\n");
            next;
         }
         $GENE = uc($GENE);
         DIELOG($outLog, "\n\ndied at processing $LCY$footPeak_logFile$N: can't parse index file def gene lqines\n\n$line\n\n") if not defined $STRAND;
         %{$coor{$GENE}} = ("CHR" => $CHR, "BEG" => $BEG, "END" => $END, "VAL" => $VAL);
         $coor{$GENE}{STRAND} = $STRAND eq "+" ? "Pos" : $STRAND eq "-" ? "Neg" : $STRAND =~ /^(Pos|Neg|Unk)$/ ? $STRAND : "Unk";
      }
      elsif ($line =~ /^-t thrshld\s+:/) {
         ($thres) = $line =~ /^-t thrshld\s+:\s+(\-?\d+\.?\d*)$/;
         $thres = "0." . $thres if $thres > 1;
      }
      elsif ($line =~ /^-w window\s+:/) {
         ($window) = $line =~ /^-w window\s+:\s+(\-?\d+\.?\d*)$/;
      }
   }
	return (\%coor, $outLog);
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

__END__
	else 
	if ($number =~ /^\-?0\.[0]*[1-9]\d*e?\-?\d*$/) {
		return($mod*abs($number));
		$number = abs($number);
		print "$number\n";
		my ($zero) = $number =~ /^0\.(0*)[1-9]\d?\d?\d*\d?\-?\d*$/;
		$zero = defined $zero ? length($zero)+1 : 1; #0.00123 = zero is 2 but times 10**5 to 123 divide by 10**2 = 1.23
		my ($scient) = $number =~ /e(\-?\d+)$/ if $number =~ /e\-?\d+$/;
		($number) = $number =~ /^(.+)e.+$/ if defined $scient;
		print "$number\n";
		$scient = 0 if not defined $scient;
		print "number=$number scient = $scient\n";
		$scient = $scient < 0 ? -1 * $zero + $scient : $zero + $scient;
		print "scient = $scient\n";
		print "$mod*abs(int($number * 10**($zero+2)+0.5)/10**2) . e$scient\n";
		return($mod*abs(int($number * 10**($zero+2)+0.5)/10**2) . "e$scient");
	}
	else {
		return $number;	
	}

