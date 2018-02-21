#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_s $opt_i $opt_g $opt_f $opt_S $opt_c $opt_C $opt_o);
getopts("s:i:g:f:S:cCo:");
BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}

use myFootLib;
my ($dupeFolder) = $opt_f;
my ($outDir) = $opt_o;
die "\nusage: $YW$0$N $CY-f [folder of -n footLop.pl]$N -n [folder of -n footPeak.pl] $LGN-o$N [output dir]\n\n" unless ex([$opt_s,$opt_S,$opt_i,$opt_g]) == 1 or ex($opt_f) == 1;
die "\nplease define output (-o)\n" if not defined $opt_o;
makedir($outDir);
my $logFile = "$dupeFolder/logFile.txt";

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
my $debugFile = "$outDir/dupefoot_debug.txt";
#my $outID = "$outDir/$samName_dupefoot.temp";
my %refz;
#open (my $outIDID, ">", "$outID") or die "Cannot write to $outID: $!\n";
open (my $outdebug, ">", "$debugFile") or die "Cannot write to $debugFile: $!\n";
my ($total_read) = `awk '\$2 == 0|| \$2 == 16 {print}' $samFile | wc -l` =~ /^(\d+)$/;
$linecount = 0;
my %read;
open (my $in1, "<", $samFile) or die "Cannot read from $samFile: $!\n";
	open (my $out1, ">", "$outDir/$samName.out");
print $out1 ">REF\n";
foreach my $chr (sort keys %refs) {
	print $out1 "$chr";
	for (my $i = 0; $i < @{$refs{$chr}}; $i+= 100) {
		my $chunk = join("", @{$refs{$chr}}[$i..$i+99]);
		my ($AT) = $chunk =~ tr/AaTt/AaTt/; $AT = 0 if not defined $AT;
		print $out1 "\t$AT";
	}
	print $out1 "\n>ENDREF\n";
}
while (my $line = <$in1>) {
	chomp($line);
	my @arr = split("\t", $line);
	next if @arr < 6;
	$linecount ++;
#	my $totalz = (keys %{$read{CALM3}});
	print date() . "\t$0: Parsed $LGN$linecount$N / $LCY$total_read$N\n" if $linecount % 50 == 0;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted, @others) = @arr;
	$chr = uc($chr);
#	next if $chr ne "CALM3";
	my $others = join("\t", @others); $others = @others == 0 ? "" : "\t$others";
	my @ref1 = defined $refs{$chr} ? @{$refs{$chr}} : die "Can't find gene $chr in $seqFile!\n";
	my $glen = @ref1;
	my ($seq2, $poz, $seqborder0, $seqborder1, $dist) = parse_samFile($line, $refs{$chr}, $out1);
}
print "Output: $outDir/$samName.out\n";
exit;

sub compare {
   my ($ref, $val1, $val2, $isPrint) = @_;
   my @ref = @{$ref};
   my ($border1L, $border1R, $valz1) = @{$val1}; my @val1 = @{$valz1};
   my ($border2L, $border2R, $valz2) = @{$val2}; my @val2 = @{$valz2};
	my $borderL = $border1L < $border2L ? $border2L : $border1L;
	my $borderR = $border1R < $border2R ? $border1R : $border2R;
	my $reflen = @ref;
	my $val1len = @val1;
	my $val2len = @val2;
	my ($print1, $print2, $print3) = ("","","");
   my ($border1, $border2) = (0,0);
   my $diff = 0; my $total = 0; my $same = 0; 
	my ($CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff) = (0,0,0,0,0,0,0,0);
   for (my $i = 0; $i < @val2; $i++) {#$borderL; $i < $borderR; $i++) {
#      if ($border1 == 0 and $val1[$i] ne "-") {
#         $border1 = 1;
#      }
#      if ($border2 == 0 and $val2[$i] ne "-") {
#         $border2 = 1;
#      }
#      if ($border1 == 1 and $val1[$i] eq "-" and join("", @val1[$i..(@val1-1)] ) =~ /^\-+$/) {
#         $border1 = 2;
#      }
#      if ($border2 == 1 and $val2[$i] eq "-" and join("", @val2[$i..(@val2-1)] ) =~ /^\-+$/) {
#         $border2 = 2;
#      }
		if ($i < $borderL or $i >= $borderR) {
			$print1 .= $val1[$i];
			$print2 .= $val2[$i];
			$print3 .= " ";
		}
      elsif ($i >= $borderL and $i < $borderR) {
			$print1 .= $val1[$i];
			$print2 .= $val2[$i];
			$print3 .= $val1[$i] eq $val2[$i] ? 1 : " ";
#=comment
			$diff ++ if $val1[$i] ne $val2[$i];
			$total ++;
			$same += $val1[$i] eq $val2[$i] ? 1 : 0;
			$CTsame += ($ref[$i] eq "C" and ($val1[$i] eq "T" or $val2[$i] eq "T") and $val1[$i] eq $val2[$i]) ? 1 : 0;
			$CTdiff += ($ref[$i] eq "C" and ($val1[$i] eq "T" or $val2[$i] eq "T") and $val1[$i] ne $val2[$i]) ? 1 : 0;
			$GAsame += ($ref[$i] eq "G" and ($val1[$i] eq "A" or $val2[$i] eq "A") and $val1[$i] eq $val2[$i]) ? 1 : 0;
			$GAdiff += ($ref[$i] eq "G" and ($val1[$i] eq "A" or $val2[$i] eq "A") and $val1[$i] ne $val2[$i]) ? 1 : 0;
			$CCsame += ($ref[$i] eq "C" and ($val1[$i] ne "T" and $val2[$i] ne "T" and ($val1[$i] eq "C" or $val2[$i] eq "C")) and $val1[$i] eq $val2[$i]) ? 1 : 0;
			$CCdiff += ($ref[$i] eq "C" and ($val1[$i] ne "T" and $val2[$i] ne "T" and ($val1[$i] eq "C" or $val2[$i] eq "C")) and $val1[$i] ne $val2[$i]) ? 1 : 0;
			$GGsame += ($ref[$i] eq "G" and ($val1[$i] ne "A" and $val2[$i] ne "A" and ($val1[$i] eq "G" or $val2[$i] eq "G")) and $val1[$i] eq $val2[$i]) ? 1 : 0;
			$GGdiff += ($ref[$i] eq "G" and ($val1[$i] ne "A" and $val2[$i] ne "A" and ($val1[$i] eq "G" or $val2[$i] eq "G")) and $val1[$i] ne $val2[$i]) ? 1 : 0;
#=cut
      }
		else {
			die "len=$val1len, $val2len, $reflen\n\nPRINT:\n1=$print1\n2=$print2\n3=$print3\n\nval2 i=$i isnt defined (val1=$val1[$i],$val1[$i+1])\n" if not defined $val2[$i]; 
			$print1 .= $val1[$i];
			$print2 .= $val2[$i];
			$print3 .= " ";
		}
		if ($i % 100 == 0) {
			#print "$print1\n$print2\n$print3\n\n" if defined $isPrint;
			$print1 = "";
			$print2 = "";
			$print3 = "";
		}
   }
	my ($L, $R) = ($borderL, $borderR);
##  return "$diff/$total";
	my $Ctot = $CCsame + $CCdiff + $CTsame + $CTdiff;
	my $Gtot = $GGsame + $GGdiff + $GAsame + $GAdiff;
##	print "CC: same=$CCsame diff=$CCdiff tot=$Ctot\n";
##	print "CT: same=$CTsame diff=$CTdiff tot=$Ctot\n";
##	print "GG: same=$GGsame diff=$GGdiff tot=$Gtot\n";
##	print "GA: same=$GAsame diff=$GAdiff tot=$Gtot\n";
##	my $Cdiff = $Ctot == 0 ? 0 : $CCdiff + $CTdiff;##int(($CCdiff + $CTdiff) / ($Ctot) * 1000+0.5)/10;
##	my $Gdiff = $Gtot == 0 ? 0 : $GGdiff + $GAdiff;##int(($GGdiff + $GAdiff) / ($Gtot) * 1000+0.5)/10;
##	print "SAME=$LGN$same$N, DIFF=$LCY$diff$N, total=$YW$total$N, Cdiff=$LGN$Cdiff$N, Ctot=$YW$Ctot$N, Gdiff=$LGN$Gdiff$N, Gtot=$YW$Gtot$N\n";
##	$Cdiff = $Ctot == 0 ? 0 : int(($CCdiff + $CTdiff) / ($Ctot) * 1000+0.5)/10;
##	$Gdiff = $Gtot == 0 ? 0 : int(($GGdiff + $GAdiff) / ($Gtot) * 1000+0.5)/10;
##   $diff = $total == 0 ? 0 : int($diff/$total*1000+0.5)/10;
##	print "len=$val1len, $val2len, $reflen\n\n";
#	my ($Ctot, $Gtot) = (0,0);
#	($diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff) = (0,0,0,0,0,0,0,0,0,0,0);
#	my @arr = (\@ref, \@val1, \@val2, $L, $R, $diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff);
#	($diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff, $Ctot, $Gtot) = compare2(@arr);##\@ref, \@val1, \@val2, $L, $R, \@arr);

##	$Ctot = $CCsame + $CCdiff + $CTsame + $CTdiff;
##	$Gtot = $GGsame + $GGdiff + $GAsame + $GAdiff;
#	print "CC: same=$CCsame diff=$CCdiff tot=$Ctot\n";
#	print "CT: same=$CTsame diff=$CTdiff tot=$Ctot\n";
#	print "GG: same=$GGsame diff=$GGdiff tot=$Gtot\n";
#	print "GA: same=$GAsame diff=$GAdiff tot=$Gtot\n";
	my $Cdiff = $Ctot == 0 ? 0 : $CCdiff + $CTdiff;#int(($CCdiff + $CTdiff) / ($Ctot) * 1000+0.5)/10;
	my $Gdiff = $Gtot == 0 ? 0 : $GGdiff + $GAdiff;#int(($GGdiff + $GAdiff) / ($Gtot) * 1000+0.5)/10;
#	print "SAME=$LGN$same$N, DIFF=$LCY$diff$N, total=$YW$total$N, Cdiff=$LGN$Cdiff$N, Ctot=$YW$Ctot$N, Gdiff=$LGN$Gdiff$N, Gtot=$YW$Gtot$N\n";
	$Cdiff = $Ctot == 0 ? 0 : int(($CCdiff + $CTdiff) / ($Ctot) * 1000+0.5)/10;
	$Gdiff = $Gtot == 0 ? 0 : int(($GGdiff + $GAdiff) / ($Gtot) * 1000+0.5)/10;
   $diff = $total == 0 ? 0 : int($diff/$total*1000+0.5)/10;
#	print "len=$val1len, $val2len, $reflen\n\n";
   return ($diff, $Cdiff, $Gdiff);
#  return int($diff/$total*1000+0.5)/10;
}

sub compare2 {
	my ($ref, $val1, $val2, $L, $R, $diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff) = @_;#[5..(@_-1)];
	#my ($diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff) = (0,0,0,0,0,0,0,0,0,0,0);
#	print "$L->$R\n";
	my $step = 25;
	my ($Ctot, $Gtot);
	for (my $i = $L; $i < $R; $i+= $step) {
		my $refz1 = join("", @{$ref}[$i..($i+$step-1)]);
		my $valz1 = join("", @{$val1}[$i..($i+$step-1)]);
		my $valz2 = join("", @{$val2}[$i..($i+$step-1)]);
		my $max = $i >= $R - $step ? $R : $i+$step;
		($Ctot) += $refz1 =~ tr/C/C/;
		($Gtot) += $refz1 =~ tr/G/G/;
#		print "$i->$max\n$valz1\n$valz2\n";
		if ($valz1 eq $valz2) {
			$total += $max - $i;#$step;
			$same += $max - $i;#$step;
			next;
#			print "\n" and next;
		}
#		print "$LRD WRONG$N\n" if $valz1 ne $valz2;
		next if $valz1 eq $valz2;
		for (my $j = $i; $j < $max; $j++) {
			$diff ++ if $val1->[$j] ne $val2->[$j];
			$total ++;
			$same += $val1->[$j] eq $val2->[$j] ? 1 : 0;
			$CTsame += ($ref->[$j] eq "C" and ($val1->[$j] eq "T" or $val2->[$j] eq "T") and $val1->[$j] eq $val2->[$j]) ? 1 : 0;
			$CTdiff += ($ref->[$j] eq "C" and ($val1->[$j] eq "T" or $val2->[$j] eq "T") and $val1->[$j] ne $val2->[$j]) ? 1 : 0;
			$GAsame += ($ref->[$j] eq "G" and ($val1->[$j] eq "A" or $val2->[$j] eq "A") and $val1->[$j] eq $val2->[$j]) ? 1 : 0;
			$GAdiff += ($ref->[$j] eq "G" and ($val1->[$j] eq "A" or $val2->[$j] eq "A") and $val1->[$j] ne $val2->[$j]) ? 1 : 0;
			$CCsame += ($ref->[$j] eq "C" and ($val1->[$j] ne "T" and $val2->[$j] ne "T" and ($val1->[$j] eq "C" or $val2->[$j] eq "C")) and $val1->[$j] eq $val2->[$j]) ? 1 : 0;
			$CCdiff += ($ref->[$j] eq "C" and ($val1->[$j] ne "T" and $val2->[$j] ne "T" and ($val1->[$j] eq "C" or $val2->[$j] eq "C")) and $val1->[$j] ne $val2->[$j]) ? 1 : 0;
			$GGsame += ($ref->[$j] eq "G" and ($val1->[$j] ne "A" and $val2->[$j] ne "A" and ($val1->[$j] eq "G" or $val2->[$j] eq "G")) and $val1->[$j] eq $val2->[$j]) ? 1 : 0;
			$GGdiff += ($ref->[$j] eq "G" and ($val1->[$j] ne "A" and $val2->[$j] ne "A" and ($val1->[$j] eq "G" or $val2->[$j] eq "G")) and $val1->[$j] ne $val2->[$j]) ? 1 : 0;
		}
	}
	return($diff, $same, $total, $CCsame, $CCdiff, $CTsame, $CTdiff, $GGsame, $GGdiff, $GAsame, $GAdiff, $Ctot, $Gtot);# = (0,0,0,0,0,0,0,0,0,0,0);
}


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
	print "Doing $seqFile\n";
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
	my $log = ""; my $log2 = "";
	my ($samFile, $seqFile, $genez, $footLoopOutDir);
	if (-e $logFile) {
		my @line = `cat $logFile`;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i]; chomp($line);
			if ($line =~ /Run script[ \t]*:[ \t]*/) {
				my ($currline) = $line =~ /Run script[ \t]*:[ \t]*(.+)$/;
				my @runscript = split(" -", $currline);
				for (my $j = 0; $j < @runscript; $j++) {
					if ($runscript[$j] =~ /^n /) {
						($footLoopOutDir) = $runscript[$j] =~ /^n (.+)$/; 
						$footLoopOutDir =~ s/\/+/\//g;
					}
				}
				my $dupeFolderTemp = getFullpath("./" . $dupeFolder);
				$log .= "Cannot parse footLoop OutDir from run script -n (line=\n\n$line\n\n)\n" and die if not defined $footLoopOutDir;
				$log .= "\n$YW-------------------->$N\n\nfootLoop.pl logFile $LGN$logFile$N parsed\n";
				$log .= "-n footLoopOutDir = $footLoopOutDir\n-n currentOutDir = $dupeFolderTemp\n\n";

			}
			if ($line =~ /^\!samFile=/) {
				($samFile) = $line =~ /^\!samFile=(.+)$/;
				$log .= "Cannot get sam File name from line=\n\n$line\n\n" and die if not defined $samFile;
				$samFile = find_file($samFile, $footLoopOutDir) if not -e $samFile;
				$log .= "Cannot get sam File name from line=\n\n$line\n\n" and die if not -e $samFile;
				$log .= "!samFile = $samFile\n";
			}
			if ($line =~ /^\!seqFile=/) {
				($seqFile) = $line =~ /^\!seqFile=(.+)$/;
				$log .= "Cannot get seq File name from line=\n\n$line\n\n" and die if not defined $seqFile;
				$seqFile = find_file($seqFile, $footLoopOutDir) if not -e $seqFile;
				$log .= "Cannot get seq File name from line=\n\n$line\n\n" and die if not -e $seqFile;
				$log .= "!seqFile = $seqFile\n";
			}
#			if ($line =~ /^!seqFile=/) {
#				($seqFile) = $line =~ /^!seqFile=(.+)$/;
#			}
			if ($line =~ /gene=.+length=/) {
				my ($gene, $length) = $line =~ /^.+gene=(.+) length=(\d+)$/;
				die if not defined $gene or not defined $length;
				$genez->{$gene} = $length;
				$log2 .= "!gene $gene = $length bp\n";
			}
#			if ($line =~ /2. Parsing in sequence for genes from sequence file/) {
#				for (my $j = $i+1; $j < @line; $j++) {
#					#$log .= "\t$line\n";
#					last if $line =~ /SUCCESS.+Sequence has been parsed from fasta file/;
#					}
#				}
#			}
			last if $line =~ /footLoop_2_sam_to_peak.pl/;
		}
	}
	print "$log\n$log2\n\n$YW<------------------$N\n\n";
	return($samFile, $seqFile, $genez);
}

sub find {
	my ($folder, $query) = @_;
	my @list = `ls -R $folder`;
	my $dir;
	my @res;
	for (my $i = 0; $i < @list; $i++) {
		chomp($list[$i]);
		if ($list[$i] =~ /:$/) {
			$dir = $list[$i];
			$dir =~ s/:$//;
		}
		else {
			my $fileName = $list[$i];
			my $file = "$dir/$list[$i]";
			if ($file =~ /$query/ or $query =~ /$file/) {
				push(@res, $file);
			}
		}
	}
	return \@res;

}

sub find_file {
	my ($samFile, $footLoopOutDir) = @_;
	if (not -e $samFile) {
		$samFile =~ s/\/+/\//g;
		if (defined $footLoopOutDir) {
			my ($samFileName) = $samFile =~ /^$footLoopOutDir(.+)$/;
			$samFile = $dupeFolder . "/" . $samFileName if defined $samFileName;
		}
		if (not -e $samFile) {
			my ($samFileName) = getFilename($samFile, "full");
			my $res = find($dupeFolder, $samFileName);
			$samFileName =~ s/\.\w+$// if $samFileName =~ /\..+\..+$/ and not defined $res;
			$res = find($dupeFolder, $samFileName) if not defined $res;
			$samFile = $res->[0] if defined $res and @{$res} == 1;
			if (defined $res and @{$res} > 1) {
				for (my $j = 0; $j < @{$res}; $j++) {
					print "$j. $res->[$j]\n";
				}
				print "Please choose:\n";
				my $choose = <STDIN>;
				while ($choose !~ /^\d+$/ or ($choose =~ /^\d+$/ and $choose >= @{$res})) {
					print "Please choose: again from the list above\n";
					$choose = <STDIN>;
				}
				$samFile = $res->[$choose];
			}
			$samFile = $res->[0] if defined $res and @{$res} == 1;
		}
	}
	print "Died at $0 parse_logFile because file ($LCY$samFile$N) does not exist in -n $LGN$dupeFolder$N\n\n" and die if not defined $samFile or (defined $samFile and not -e $samFile);
	$samFile = getFullpath($samFile);
	return $samFile;
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
	my ($line, $refs, $out) = @_;
	my %dist;
	my ($read, $strand, $chr, $pos, $mapq, $cigar, $junk1, $junk2, $junk3, $seqs, $qual, $junk4, $junk5, $converted) = split("\t", $line);
	#print "POS=$pos\n";
	my @seq1 = split("", $seqs);
	my @ref1 = @{$refs};
	my (@ref2, @seq2);
	my ($refpos, $seqpos) = ($pos-1, 0);
	my $seqborder1 = $pos-1;
	#@seq2 = ("-")x($pos-1);
	#@ref2 = @ref1[0..$pos-2];
	##print join("", @ref2) . "\n";
	##print join("", @seq2) . "\n\n";
	my ($num, $alp, $lengthseq) = parse_cigar($cigar); die if not defined $num;
	for (my $i = 0; $i < @{$num}; $i++) {
		if ($alp->[$i] eq "I") {
#			$dist{$refpos . ".5"} = "I";#.$i.$num->[$i]";
			$dist{$refpos . ".5"} = "I.$num->[$i]" if $ref1[$refpos] !~ /[CG]/;
			$seqpos += $num->[$i];
		}
		elsif ($alp->[$i] eq "D") {
			#@ref2 = (@ref2, @ref1[$refpos..($refpos+$num->[$i]-1)]);
			#@seq2 = (@seq2, ("-")x$num->[$i]);
			$seqborder1 += $num->[$i];
#			$dist{$refpos} = "D";#"D.$i.$num->[$i]";
			$dist{$refpos} = "D.$num->[$i]" if $ref1[$refpos] !~ /[CG]/;
#			for (my $k = $refpos; $k < $refpos+$num->[$i]; $k++) {
#				$dist{$k} = "D.$i.$num->[$i]";
#			}
			$refpos += $num->[$i];
		}
		else {
			my ($begseq, $endseq) = ($seqpos, $seqpos+$num->[$i]);
			#@seq2 = (@seq2, @seq1[$begseq..$endseq-1]);
			$seqborder1 += ($endseq-$begseq);
			$endseq = $endseq >= @seq1 ? @seq1-1 : $endseq;
			#@ref2 = (@ref2, @ref1[$refpos..($refpos+$num->[$i]-1)]); 
##			print "begref=$refpos, begseq=$begseq, begrefnuc=\n" . join("", @ref1[$refpos..$refpos+10]) . "\n" . join("", @seq1[$begseq..$begseq+10]) . "\n" if $refpos <= 100;
			if (join("", @seq1[$begseq..$endseq-1]) ne join("", @ref1[$refpos..($refpos+$num->[$i]-1)])) {
				my ($k, $l) = ($begseq, $refpos);
				while ($k < $endseq) {
					if ($ref1[$l] ne $seq1[$k]) {
						my $seqz = ($ref1[$l] eq "C" and $seq1[$k] eq "T") ? "X" : ($ref1[$l] eq "G" and $seq1[$k] eq "A") ? "Y" : "$seq1[$k]";
						if ($seqz =~ /[XY]/ or $ref1[$l] !~ /[CG]/) {
							$dist{$l} = "$seqz";
						}
						#$dist{$l} = ($ref1[$l] eq "C" and $seq1[$k] eq "T") ? "X" : ($ref1[$l] eq "G" and $seq1[$k] eq "A") ? "Y" : "$seq1[$k]";
					}
					$k ++; $l ++;
				}
			}
			($refpos, $seqpos) = ($refpos + $num->[$i], $seqpos + $num->[$i]);
		}
	}
	my $seqborder0 = $pos - 1;
	#my $seqborder1 = @seq2;
	#@ref2 = (@ref2, @ref1[@ref2..@ref1-1]);
	#@seq2 = (@seq2, ("-")x(@ref2-@seq2));
	print $out1 "$read,$chr,$mapq,$seqborder0,$seqborder1";
	foreach my $pos (sort {$a <=> $b} keys %dist) {
		print $out1 ",$pos=$dist{$pos}";
	}
	print $out1 "\n";
#	print join("", @ref1) . "\n";
#	print join("", @seq1) . "\n";
##	if ($seqborder1 >= 2900) {
#	print "REF=," . scalar(@ref1) . ",2=" . scalar(@ref2) . "\n" . join("", @ref2) . "\nSEQ=" . scalar(@seq2) . "\n" . join("", @seq2) . "\n\n";
#	for (my $i = 0; $i < @ref2; $i++) {
#		print " " if not defined $dist{$i};
#		print "$dist{$i}" if defined $dist{$i};
#	}
#	print "\n";
#	for (my $i = 0; $i < @ref2; $i++) {
#		print " " if not defined $dist{$i . ".5"};
#		print $dist{$i . ".5"} if defined $dist{$i . ".5"};
#	}
##	}
##	compare(\@ref2, [$seqborder0, $seqborder1, \@ref2], [$seqborder0, $seqborder1, \@seq2], 1);#$seqborder0, $seqborder1, 1);
#	die;
	#return(\@seq2, 1, $seqborder0, $seqborder1, \%dist);#$seqborder0, $seqborder1);
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
	my (@ref2, @seq2);
	for (my $i = 0; $i < @ref; $i++) {
      if ($i < $seqborder0 or $i >= $seqborder1) {
			push(@ref2, $ref[$i]);
			push(@seq2, $seq[$i]);
      }
		else {
			push(@ref2, $ref[$i]) if $ref[$i] ne "-";
			push(@seq2, $seq[$i]) if $ref[$i] ne "-";
		}
	}
	#compare(\@ref2, \@ref2, \@seq2, 1);
	return(\@ref2, \@seq2, \%pos, $seqborder0, $seqborder1);
}
