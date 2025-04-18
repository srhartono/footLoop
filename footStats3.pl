#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_d $opt_n $opt_G $opt_c $opt_F);
getopts("vdn:G:cF");
use mitochy;

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
	print "\n- Pushed $libPath into perl lib path INC\n";

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
}

use myFootLib;
use FAlite;

my $homedir = $ENV{"HOME"};

my ($footLoop_script_folder, $version, $md5script) = check_software();
my ($version_small) = $version =~ /^([vV]\d+\.\d+)[a-zA-Z_]*.*$/;
$version_small = $version if not defined $version_small;

#my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
#my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
#my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
#my $version = join("", @version);
#if (defined $opt_v) {
#   print "$version\n";
#   exit;
#}
#my ($version_small) = "vUNKNOWN";
#foreach my $versionz (@version[0..@version-1]) {
#   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
#}

my $date = getDate();
my $uuid = getuuid();
my $numThread = 1;
my ($footPeakFolder, $outLog) = check_sanity();
my $genewant = $opt_G;
my $opts = footStats_parse_footPeak_logFile($footPeakFolder, "$footPeakFolder/footPeak_logFile.txt", $outLog);
my ($footLoopFolder) = $opts->{footLoop2}{n};
my $geneIndexFile = $opts->{footLoop2}{i};
my $coor = parse_geneIndexFile($geneIndexFile, $outLog);
DIELOG($outLog, date() . " ERROR: $footLoopFolder footloop folder doesn't exist!\n") if not -d $footLoopFolder;
my %data;
my ($labelFile) = "$footPeakFolder/.LABEL";
my @nuc = qw(0 1 2 3 4 5 6 7 8 9);
my @color = ($GR, $N, $N, $N, $YW, $YW, $LGN, $LGN, $LRD, $LRD); 

#for (my $i = 0; $i < @peakFiles; $i++) {
#	my $peakFile = $peakFiles[$i];
#	print "$i $peakFile\n";
#	my $file = $peakFile;
#	print "file=$file\n";
#	open (my $in1, "<", $file) or LOG($outLog, "Failed to read from $file: $!\n") and exit 1;
#	while (my $line = <$in1>) {
#		chomp($line);
#		my ($read, @val) = split("\t", $line);
#		print "HERE\n$read\n";
#		last;
#	}
#	close $in1;
#	last if $i > 10;
#}
my $wclFile = "$footPeakFolder/99_FOOTSTATS/0_wcl.tsv";
if (not -e $wclFile or defined $opt_F) {
	my $cmd = "footPeak_wclthese.pl -n $footPeakFolder";
	#print "\nGetting linecount of each .out files cmd #1:\n$LCY$cmd$N\n\n";
	system($cmd) == 0 or die "Failed to run cmd: $!\n$LCY$cmd$N\n\n";
	my $cmd2 = "cat $footPeakFolder/.CALL/*.wcl > $footPeakFolder/99_FOOTSTATS/0_wcl.tsv";
	#print "\nGetting linecount of each .out files cmd #2:\n$LCY$cmd2$N\n\n";
	system($cmd2) == 0 or die "Failed to run cmd: $!\n$LCY$cmd2$N\n\n";
}
my $wclhash = parse_wcl($wclFile);
my %wcl = %{$wclhash};
sub parse_wcl {
	my ($wclFile) = @_;
	my %wcl;
	my $linecount = 0;
	open (my $in1, "<", $wclFile) or die "Failed to read from $wclFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		$line =~ s/^\s+(\d+)/$1/;
		$line =~ s/\/+/\//g;
		my ($totalread, $file) = split(" ", $line);
		$file = "$footPeakFolder/.CALL/$file" if $file !~ /\/\.CALL\//;
		$file =~ s/\/+/\//g;
		$file =~ s/\.\//\//g;
		$file =~ s/\/+/\//g;
		#print "\n${LGN}wclFile $LCY$wclFile$N: Example line $linecount:\n$file\t$totalread\n\n" if $linecount eq 1;
		$wcl{$file} = $totalread;
	}
	close $in1;
	return(\%wcl);
}

my $peakFileCount = 0;
my $totalfile = (keys %wcl);
my %bed;
my $lastplasmid;
print "Checking peak.out Files (n=$LGN$totalfile$N)\n\n";
$peakFileCount = 0;
my $counthash;
print "${LCY}PBEH2_BC81_PLASMIDPFC8_TACT1T2_DESCLINEAR_SGRNA08_TX.m84066_240320_204128_s1/123669377/ccs$N\t${LPR}peaktype$N\tconvperc\tconv\tlength";
foreach my $nuc (@nuc[0..@nuc-1]) {
	print "\t$color[$nuc]$nuc     $N";
}
print "\n";
foreach my $peakFile (sort keys %wcl) {
	if (defined $opt_G and $peakFile !~ /$opt_G/i) {
		next;
	}
	$peakFileCount ++;
	my $totalread = $wcl{$peakFile};
	$data{$peakFile}{PEAK}{totalread} = $totalread;
	#print "$LGN$peakFileCount/$totalfile$N: totalread=$YW$totalread$N $LCY$bedFile$N\n" if -e $bedFile;next;
	print "$peakFileCount/$totalfile$N $peakFile\n" if $peakFileCount % 1000 == 0;
	print "$LGN$peakFileCount/$totalfile$N totalread=$YW$totalread$N\n" if $peakFileCount % 1000 == 0;
	my $parseName = parseName($peakFile);
	my ($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	my $readstrand = $strand;
	$label =~ s/_bismark_bt2.bam//;
	my $plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	my $genestrand = $coor->{$plasmid2}{strand}; 
	$genestrand = "Pos" if not defined $genestrand;
	$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	my $array = $parseName->{$peakFile}{array};
   my ($flag) = getFlag($peakFile, $genestrand, $readstrand, $convtype);
	#print "$LCY$plasmid2$N $YW$readstrand$N $LGN$convtype$N $LPR$flag$N\n";
	$lastplasmid = $plasmid2 if not defined $lastplasmid;

	if ($lastplasmid ne $plasmid2) {
		print "lastplasmid=$LGN$lastplasmid$N plasmid2=$LCY$plasmid2$N\n";
		printcounthash($counthash);
		undef $counthash;
	}
	$lastplasmid = $plasmid2;
	$counthash = parse_outFile($peakFile, $counthash, $flag);# if $flag eq "PEAK_C";
	#print "$lastplasmid $plasmid2 total $flag = " . (keys %{$counthash}) . "\n";
	last if $peakFileCount > 320;
}

		#print "Here\n";
if (1 == 1) {
	print "lastplasmid=$LGN$lastplasmid$N\n";
	printcounthash($counthash);
}
sub printcounthash {
	my ($counthash) = @_;
		#my %printed;
		my $totalcounthash = (keys %{$counthash});
		print "$totalcounthash\n";
		foreach my $read (sort keys %{$counthash}) {
			next if not defined $counthash->{$read};
			next if (keys %{$counthash->{$read}} < 4);
			next if not defined $counthash->{$read}{NOPK_TEMP_C};
			next if not defined $counthash->{$read}{NOPK_TEMP_RCONV_C};
		#	$printed{$flag} = 1;
#			print "$LCY$read$N";
			my $perc0  = $counthash->{$read}{NOPK_TEMP_C}{perc}{6} + $counthash->{$read}{NOPK_TEMP_C}{perc}{7} + $counthash->{$read}{NOPK_TEMP_C}{perc}{8} + $counthash->{$read}{NOPK_TEMP_C}{perc}{9};
			my $total0 = $counthash->{$read}{NOPK_TEMP_C}{perc}{4} + $counthash->{$read}{NOPK_TEMP_C}{perc}{5} + $perc0;
			my $len0 = int($counthash->{$read}{NOPK_TEMP_C}{total} * $total0/100+0.5);
			my $finalperc0 = $total0 == 0 ? 1 : int(1000*$perc0 / $total0+0.5)/10;
			my $conv0 = int($counthash->{$read}{NOPK_TEMP_C}{total} * $perc0/100+0.5);

			my $perc1  = $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{6} + $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{7} + $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{8} + $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{9};
			my $total1 = $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{4} + $counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{5} + $perc1;
			my $len1 = int($counthash->{$read}{NOPK_TEMP_RCONV_C}{total} * $total1/100+0.5);
			my $finalperc1 = $total1 == 0 ? 1 : int(1000*$perc1 / $total1+0.5)/10;
			my $conv1 = int($counthash->{$read}{NOPK_TEMP_RCONV_C}{total} * $perc1/100+0.5);

			print "$LCY$read$N\t${LPR}NOPK_TEMP_C$N\t$finalperc0\t$conv0\t$len0";
			foreach my $nuc (@nuc[0..@nuc-1]) {
				print "\t$color[$nuc]$counthash->{$read}{NOPK_TEMP_C}{perc}{$nuc}$N";
			}
			print "\n";
			print "$LCY$read$N\t${LPR}NOPK_TEMP_RCONV_C$N\t$finalperc1\t$conv1\t$len1";
			foreach my $nuc (@nuc[0..@nuc-1]) {
				print "\t$color[$nuc]$counthash->{$read}{NOPK_TEMP_RCONV_C}{perc}{$nuc}$N";
			}
			print "\n";
		#	last if (keys %printed == 16);
		}
		print "\n\n";

		foreach my $read (sort keys %{$counthash}) {
			next if not defined $counthash->{$read};
			next if (keys %{$counthash->{$read}} < 4);
			next if not defined $counthash->{$read}{NOPK_C};
			next if not defined $counthash->{$read}{NOPK_RCONV_C};
#			$printed{$flag} = 1;
#			print "$LCY$read$N" if $flag =~ /_C/;
			my $perc0  = $counthash->{$read}{NOPK_C}{perc}{6} + $counthash->{$read}{NOPK_C}{perc}{7} + $counthash->{$read}{NOPK_C}{perc}{8} + $counthash->{$read}{NOPK_C}{perc}{9};
			my $total0 = $counthash->{$read}{NOPK_C}{perc}{4} + $counthash->{$read}{NOPK_C}{perc}{5} + $perc0;
			my $len0 = int($counthash->{$read}{NOPK_C}{total} * $total0/100+0.5);
			my $finalperc0 = $total0 == 0 ? 0 : int(1000*$perc0 / $total0+0.5)/10;
			my $conv0 = int($counthash->{$read}{NOPK_C}{total} * $perc0/100+0.5);

			my $perc1  = $counthash->{$read}{NOPK_RCONV_C}{perc}{6} + $counthash->{$read}{NOPK_RCONV_C}{perc}{7} + $counthash->{$read}{NOPK_RCONV_C}{perc}{8} + $counthash->{$read}{NOPK_RCONV_C}{perc}{9};
			my $total1 = $counthash->{$read}{NOPK_RCONV_C}{perc}{4} + $counthash->{$read}{NOPK_RCONV_C}{perc}{5} + $perc1;
			my $len1 = int($counthash->{$read}{NOPK_RCONV_C}{total} * $total1/100+0.5);
			my $finalperc1 = $total1 == 0 ? 1 : int(1000*$perc1 / $total1+0.5)/10;
			my $conv1 = int($counthash->{$read}{NOPK_RCONV_C}{total} * $perc1/100+0.5);
			print "$LCY$read$N\t${LPR}NOPK_C$N\t$finalperc0\t$conv0\t$len0";
			foreach my $nuc (@nuc[0..@nuc-1]) {
				print "\t$color[$nuc]$counthash->{$read}{NOPK_C}{perc}{$nuc}$N";
			}
			print "\n";
			print "$LCY$read$N\t${LPR}NOPK_RCONV_C$N\t$finalperc1\t$conv1\t$len1";
			foreach my $nuc (@nuc[0..@nuc-1]) {
				print "\t$color[$nuc]$counthash->{$read}{NOPK_RCONV_C}{perc}{$nuc}$N";
			}
			print "\n";
#			last if (keys %printed == 16);
		}
		print "\n\n";
		undef $counthash;
}
print "DEBUG\n";
exit 0;
sub parse_outFile {
	my @flags = qw(NOPK_C NOPK_TEMP_C NOPK_RCONV_C NOPK_TEMP_RCONV_C);
	my @nuc = qw(0 1 2 3 4 5 6 7 8 9);
	my ($peakoutFile, $counthash, $flag) = @_;
	my $linecount = 0;
	open (my $in, "<", $peakoutFile) or die;
	while (my $line = <$in>) {
		chomp($line);
		my ($read, @val) = split("\t", $line);
		#foreach my $flags (@flags) {
		#	foreach my $nuc (@nuc[0..@nuc-1]) {
		#		$counthash->{$read}{$flags}{perc}{$nuc} = 0 if not defined $counthash->{$read}{$flags}{perc}{$nuc};
		#	}
		#}
		my $val = join("", @val);
		#print "$read";# if $linecount == 0;
		$linecount ++;
		$counthash->{$read}{$flag}{total} = @val;
		for (my $i = 0; $i < @val; $i++) {
			foreach my $nuc (@nuc[0..@nuc-1]) {
				$counthash->{$read}{$flag}{count}{$nuc} = 0 if not defined $counthash->{$read}{$flag}{count}{$nuc};
				$counthash->{$read}{$flag}{count}{$nuc} ++ if $val[$i] eq $nuc;
			}
			#print "\t$nuc[$i]\t$nuc[4]\t$n4\t$n42\t$counthash->{$read}{$flag}{total}\ttot=$tot\n";
			#exit 0;
		}
		foreach my $nuc (@nuc[0..@nuc-1]) {
			$counthash->{$read}{$flag}{perc}{$nuc} = $counthash->{$read}{$flag}{count}{$nuc} / $counthash->{$read}{$flag}{total};
			my $mult = $counthash->{$read}{$flag}{perc}{$nuc} < 0.01 ? 100000 : 1000;
			$counthash->{$read}{$flag}{perc}{$nuc} = int($mult*$counthash->{$read}{$flag}{perc}{$nuc}+0.5)/($mult/100);
			
		}
		foreach my $nuc (@nuc[0..@nuc-1]) {
			my @color = ($GR, $N, $N, $N, $YW, $YW, $LGN, $LGN, $LRD, $LRD); 
#			print "\t$color[$nuc]$nuc$N=
			#print "\t$color[$nuc]$counthash->{$read}{$flag}{count}{$nuc}$N";
			#print "\t$color[$nuc]$counthash->{$read}{$flag}{perc}{$nuc}$N";
		}
		#print "\n";
		last if $linecount > 10;
	}
	close $in;
	return($counthash);
}
print "Checking bedFiles (n=$LGN$totalfile$N)\n\n";
foreach my $peakFile (sort keys %wcl) {
	$data{$peakFile}{PEAKLEN}{MEAN} = 0;
	$data{$peakFile}{PEAKLEN}{MEDIAN} = 0;
	$data{$peakFile}{PEAKLEN}{SD} = 0;
	$data{$peakFile}{PEAK}{total} = 0;
	for (my $i = 0; $i <= 100; $i+=25) {
		$data{$peakFile}{PEAK}{$i} = 0;
		$data{$peakFile}{PEAKPERC}{$i} = 0;
		$data{$peakFile}{PEAKATLEAST}{$i} = 0;
		$data{$peakFile}{PEAKATLEASTPERC}{$i} = 0;
	}
	$data{$peakFile}{PEAK}{totalpeak} = 0;
	$data{$peakFile}{PEAK}{totalread} = 0;
	$peakFileCount ++;
	my ($folder1, $peakFileName1) = getFilename($peakFile, "folderfull");
	
	my $bedFile = $peakFile;
	
	$bedFile =~ s/\/\.CALL\//\/PEAKS_LOCAL\//;
	$bedFile =~ s/.out$/.local.bed/;
	my $totalread = $wcl{$peakFile};
	$data{$peakFile}{PEAK}{totalread} = $totalread;
	#print "$LGN$peakFileCount/$totalfile$N: totalread=$YW$totalread$N $LCY$bedFile$N\n" if -e $bedFile;next;
	print "$peakFileCount/$totalfile $LCY$bedFile$N\n" if $peakFileCount % 1000 == 0;
	next if not -e $bedFile;
	print "$LGN$peakFileCount/$totalfile$N totalread=$YW$totalread$N bedFile=$LCY$bedFile$N\n" if $peakFileCount % 1000 == 0;
	my $parseName = parseName($bedFile);
	my ($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	my $readstrand = $strand;
	$label =~ s/_bismark_bt2.bam//;
	my $plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	my $genestrand = $coor->{$plasmid2}{strand}; 
	$genestrand = "Pos" if not defined $genestrand;
	$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	my $array = $parseName->{$peakFile}{array};
   my ($flag) = getFlag($bedFile, $genestrand, $readstrand, $convtype);
	open (my $in1, "<", $bedFile) or die;
	my %temp;
	@{$temp{i}} = ();
	@{$temp{read}} = ();
	@{$temp{beg}} = ();
	@{$temp{end}} = ();
	@{$temp{len}} = ();
	my $totalpeak = 0;
	my $linecount = 0;
	my %lenmax;
	my %totalread;
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my ($read, $beg, $end) = split("\t", $line);
		my $len = $end - $beg;
		push(@{$temp{i}}, $linecount);
		push(@{$temp{read}}, $read);
		push(@{$temp{beg}}, $beg);
		push(@{$temp{end}}, $end);
		push(@{$temp{len}}, $len);
		#$temp[$linecount]{read} = $read;
		#$temp[$linecount]{beg} = $beg;
		#$temp[$linecount]{end} = $end;
		#$temp[$linecount]{len} = $len;
		$totalread{$read} = 1;
		$totalpeak ++;
		$lenmax{$read} = $len if not defined $lenmax{$read};
		$lenmax{$read} = $len if $len > $lenmax{$read};
	}
	close $in1;
	$totalpeak = (keys %totalread);
	$data{$peakFile}{PEAK}{totalpeak} = $totalpeak;
	#my $totalread = (keys %totalread);
	#my $sd = @{$temp{len}} == 0 ? 0 : int(1*mitochy::sd(\@{$temp{len}})+0.5)/1;
	#my $mean = @{$temp{len}} == 0 ? 0 : int(1*mitochy::mean(\@{$temp{len}})+0.5)/1;
	#my $median = @{$temp{len}} == 0 ? 0 : int(1*mitochy::median(\@{$temp{len}})+0.5)/1;
	my @lenmax;
	foreach my $read (sort keys %lenmax) {
		push(@lenmax, $lenmax{$read});
	}
	my $sd = @lenmax == 0 ? 0 : int(1*mitochy::sd(\@lenmax)+0.5)/1;
	my $mean = @lenmax == 0 ? 0 : int(1*mitochy::mean(\@lenmax)+0.5)/1;
	my $median = @lenmax == 0 ? 0 : int(1*mitochy::median(\@lenmax)+0.5)/1;
	
	$data{$peakFile}{PEAKLEN}{MEAN} = $mean;
	$data{$peakFile}{PEAKLEN}{SD} = $sd;
	$data{$peakFile}{PEAKLEN}{MEDIAN} = $median;
	my ($temphash2, $totaltemp2, %temp2);
	my @ops;my @ops2;
	my ($colname, $op, $num2);
	
	for (my $i = 0; $i <= 100; $i+=25) {
		@ops = ();
		my $i0 = $i;
		my $i1 = $i == 100 ? 999999 : $i+25;
		push(@ops, ["len","gte",$i0]);
		push(@ops, ["len","lt",$i1]);
		($temphash2, $totaltemp2) = check0(\%temp, \@ops);
		%temp2 = defined $temphash2 ? %{$temphash2} : ();
		my %total2;
		$totaltemp2 = 0;
		if (defined $temp2{i}) {
			for (my $k = 0; $k < @{$temp2{i}};$k++ ){
				my $read = $temp2{read}[$k];
				$total2{$read} ++;
			}
			$totaltemp2 = (keys %total2);
		} else {$totaltemp2 = 0;}
		my $totaltemp2perc = $totalread == 0 ? 0 : int(1000*$totaltemp2 / $totalread+0.5)/10;
		print "[$i0-$i1 " . join(" ", @{$ops[0]}) . " AND " . join(" ", @{$ops[1]}) . "]: $LGN$totaltemp2$N/$YW$totalread$N ($LGN$totaltemp2perc$N %)\n" if defined $opt_d;
		$data{$peakFile}{PEAK}{$i} = $totaltemp2;
		$data{$peakFile}{PEAKPERC}{$i} = $totaltemp2perc;
	}

	for (my $i = 0; $i <= 100; $i+=25) {
		@ops = ();
		my $i0 = $i;
		my $i1 = $i == 100 ? 999999 : $i+25;
		push(@ops, ["len","gte",$i0]);
		($temphash2,$totaltemp2) = check0(\%temp, \@ops);
		%temp2 = defined $temphash2 ? %{$temphash2} : ();
		my %total2;
		$totaltemp2 = 0;
		if (defined $temp2{i}) {
			for (my $k = 0; $k < @{$temp2{i}};$k++ ){
				my $read = $temp2{read}[$k];
				$total2{$read} ++;
			}
			$totaltemp2 = (keys %total2);
		} else {$totaltemp2 = 0;}
		my $totaltemp2perc = $totalread == 0 ? 0 : int(1000*$totaltemp2 / $totalread+0.5)/10;
		print "[$i0 " . join(" ", @{$ops[0]}) . "]: $LGN$totaltemp2$N/$YW$totalread$N ($LGN$totaltemp2perc$N %)\n" if defined $opt_d;
		$data{$peakFile}{PEAKATLEAST}{$i} = $totaltemp2;
		$data{$peakFile}{PEAKATLEASTPERC}{$i} = $totaltemp2perc;
	}
	
	#$data{$peakFile}{"$colname;$op;$num2"} = $totaltemp2;
	#print "]: $LGN$totaltemp2$N/$YW$totalpeak$N\n";
	#($temphash2, $totaltemp2) = check(\%temp, $colname, $op, $num2);
	#if ($totalpeak > 0) {
		#($temphash2, $totaltemp2) = check(\%temp, "len", "gte", 100);
		#%temp2 = defined $temphash2 ? %{$temphash2} : ();
		#print "[len >= 100]: $LGN$totaltemp2$N/$YW$totalpeak$N\n";
		#($temphash2, $totaltemp2) = check($temphash2, "len", "gte", 120);
		#%temp2 = defined $temphash2 ? %{$temphash2} : ();
		#print "[len >= 120]: $LGN$totaltemp2$N/$YW$totalpeak$N\n";
	#}
	#else {
	#	print "-> totalpeak=$LGN$totalpeak$N\n";
	#}
	last if $peakFileCount > 100 and defined $opt_d;
}
print "\n";
foreach my $peakFile (sort keys %data) {
	my $totalpeak = $data{$peakFile}{PEAK}{totalpeak};
	my $totalread = $data{$peakFile}{PEAK}{totalread};
	my $peak0 = $data{$peakFile}{PEAK}{0};
	my $peak25 = $data{$peakFile}{PEAK}{25};
	my $peak50 = $data{$peakFile}{PEAK}{50};
	my $peak75 = $data{$peakFile}{PEAK}{75};
	my $peak100 = $data{$peakFile}{PEAK}{100};
	my $peakperc0 = $data{$peakFile}{PEAKPERC}{0};
	my $peakperc25 = $data{$peakFile}{PEAKPERC}{25};
	my $peakperc50 = $data{$peakFile}{PEAKPERC}{50};
	my $peakperc75 = $data{$peakFile}{PEAKPERC}{75};
	my $peakperc100 = $data{$peakFile}{PEAKPERC}{100};
	my $peakatleast0 = $data{$peakFile}{PEAKATLEAST}{0};
	my $peakatleast25 = $data{$peakFile}{PEAKATLEAST}{25};
	my $peakatleast50 = $data{$peakFile}{PEAKATLEAST}{50};
	my $peakatleast75 = $data{$peakFile}{PEAKATLEAST}{75};
	my $peakatleast100 = $data{$peakFile}{PEAKATLEAST}{100};
	my $peakatleastperc0 = $data{$peakFile}{PEAKATLEASTPERC}{0};
	my $peakatleastperc25 = $data{$peakFile}{PEAKATLEASTPERC}{25};
	my $peakatleastperc50 = $data{$peakFile}{PEAKATLEASTPERC}{50};
	my $peakatleastperc75 = $data{$peakFile}{PEAKATLEASTPERC}{75};
	my $peakatleastperc100 = $data{$peakFile}{PEAKATLEASTPERC}{100};
#	my $totalpeak = $data{$peakFile}{PEAK}{total};
	my $mean = $data{$peakFile}{PEAKLEN}{MEAN};
	my $median = $data{$peakFile}{PEAKLEN}{MEDIAN};
	my $sd = $data{$peakFile}{PEAKLEN}{SD};
	$data{$peakFile}{PEAKPRINT} = "\t$totalread\t$totalpeak\t$mean\t$median\t$sd";
	$data{$peakFile}{PEAKPRINT} .= "\t$LCY$peak0\t$peak25\t$peak50\t$peak75\t$peak100";
	$data{$peakFile}{PEAKPRINT} .= "\t$LGN$peakatleast0\t$peakatleast25\t$peakatleast50\t$peakatleast75\t$peakatleast100";
	$data{$peakFile}{PEAKPRINT} .= "\t$LCY$peakperc0\t$peakperc25\t$peakperc50\t$peakperc75\t$peakperc100";
	$data{$peakFile}{PEAKPRINT} .= "\t$LGN$peakatleastperc0\t$peakatleastperc25\t$peakatleastperc50\t$peakatleastperc75\t$peakatleastperc100$N";
	print "$peakFile$data{$peakFile}{PEAKPRINT}\n" if defined $opt_d;
}
#	my @op;
#	for (my $i = 100; $i < 200; $i+= 100) {
#		$colname = "len";
#		$op = "gte";
#		$num2 = $i;
#		my @op2 = ($colname, $op, $num2);
#		push(@op, \@op2);
#	}
#	($temphash2, $totaltemp2) = check0(\%temp, \@op);
print "Checking .CALL/*PEAK*.out\n";
my @peakFiles = <$footPeakFolder/.CALL/*PEAK*out>;
my %final;
my ($len) = $peakFiles[0] =~ /^.+_l(\d+)/;
for (my $i = 0; $i < @peakFiles; $i++) {
	print "Done $i\n" if $i % 1000 == 0;
	my $peakFile = $peakFiles[$i];
	$peakFile =~ s/\/+/\//g;
	if ($i == 0) {
		#print "Example: \npeakFile=$LCY$peakFile$N, total=$wcl{$peakFile}\n";
	}
	($data{$peakFile}{total}) = defined($wcl{$peakFile}) ? $wcl{$peakFile} : 0;
	#($data{$peakFile}{total}) = `wc -l $peakFile` =~ /^(\d+)/ if -e $peakFile;
	my $parseName = parseName($peakFile);
	my ($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	my $readstrand = $strand;
	$label =~ s/_bismark_bt2.bam//;
	my $plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	my $genestrand = $coor->{$plasmid2}{strand}; 
	$genestrand = "Pos" if not defined $genestrand;
	$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	my $array = $data{$peakFile}{array};
   my ($flag) = getFlag($peakFile, $genestrand, $readstrand, $convtype);
	#print "$i $YW$bc $plasmid2 $plasmid $desc $LPR$genestrand $LCY$readstrand $convtype $LPR$flag$N\t$LGN$data{$peakFile}{total}$N\t$peakFile\n";

	my $nopkFile = $peakFile; $nopkFile =~ s/PEAK/NOPK/g;
	($data{$nopkFile}{total}) = defined($wcl{$nopkFile}) ? $wcl{$nopkFile} : 0;
	my $peakprint = $data{$peakFile}{PEAKPRINT};
	my $nopkprint = $data{$nopkFile}{PEAKPRINT};
	#($data{$nopkFile}{total}) = `wc -l $nopkFile` =~ /^(\d+)/ if -e $nopkFile;
	my $totalread = $data{$peakFile}{total} + $data{$nopkFile}{total};
	my $print = "$bc\t$plasmid\t$desc\t$readstrand\t$convtype\t$data{$peakFile}{total}\t$totalread\t$flag\t$peakFile\n";
	my $plasmiddesc = "PLASMID$plasmid;DESC$desc";
	#print "$print\n" if $i == 0;
	my $totalpeak = $data{$peakFile}{PEAK}{total};
	my $peakmean = $data{$peakFile}{PEAKLEN}{MEAN}; $peakmean = 0 if not defined $peakmean;
	my $nopkmean = $data{$nopkFile}{PEAKLEN}{MEAN}; $nopkmean = 0 if not defined $nopkmean;
	my $median = $data{$peakFile}{PEAKLEN}{MEDIAN};
	my $sd = $data{$peakFile}{PEAKLEN}{SD};
	my $peakatleastperc0 = $data{$peakFile}{PEAKATLEASTPERC}{0};
	my $peakatleastperc25 = $data{$peakFile}{PEAKATLEASTPERC}{25};
	my $peakatleastperc50 = $data{$peakFile}{PEAKATLEASTPERC}{50};
	my $peakatleastperc75 = $data{$peakFile}{PEAKATLEASTPERC}{75};
	my $peakatleastperc100 = $data{$peakFile}{PEAKATLEASTPERC}{100};
	my $nopkatleastperc0 = $data{$nopkFile}{PEAKATLEASTPERC}{0};
	my $nopkatleastperc25 = $data{$nopkFile}{PEAKATLEASTPERC}{25};
	my $nopkatleastperc50 = $data{$nopkFile}{PEAKATLEASTPERC}{50};
	my $nopkatleastperc75 = $data{$nopkFile}{PEAKATLEASTPERC}{75};
	my $nopkatleastperc100 = $data{$nopkFile}{PEAKATLEASTPERC}{100};
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{flag} = $flag;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak} = $data{$peakFile}{total};
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk} = $data{$nopkFile}{total};
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakperc} = $totalread == 0 ? 0 : int(1000*$data{$peakFile}{total}/$totalread+0.5)/10;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkperc} = $totalread == 0 ? 0 : int(1000*$data{$nopkFile}{total}/$totalread+0.5)/10;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread} = $totalread;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakmean} = $peakmean;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkmean} = $nopkmean;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak100} = $peakatleastperc100;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak50} = $peakatleastperc50;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak75} = $peakatleastperc75;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak25} = $peakatleastperc25;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk100} = $nopkatleastperc100;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk50} = $nopkatleastperc50;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk25} = $nopkatleastperc25;
	$final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk75} = $nopkatleastperc75;
	#$final{$bc}{$plasmiddesc}{$strand}{$convtype}{median} = $median;
	#$parseName = parseName($nopkFile);
	#($label, $gene, $strand, $window, $thres, $convtype, $bc, $plasmid, $desc, $pcb) = @{$parseName->{array}};
	#$readstrand = $strand;
	#$label =~ s/_bismark_bt2.bam//;
	#$plasmid2 = "$label\_BC$bc\_PLASMID$plasmid\_DESC$desc";
	#$genestrand = $coor->{$plasmid2}{strand}; 
	#$genestrand = "Pos" if not defined $genestrand;
	#$genestrand = $genestrand eq "+" ? "Pos" : $genestrand eq "-" ? "Neg" : $genestrand;
	#$array = $data{$nopkFile}{array};
   #($flag) = getFlag($nopkFile, $genestrand, $readstrand, $convtype);
	#print "$i\t$YW$bc\t$plasmid\t$desc\t$readstrand\t$convtype\t$LGN$data{$peakFile}{total}$N\t$LGN$totalread$N\t$LPR$flag$N\n";
	#print "$i $YW$bc $plasmid2 $plasmid $desc $LPR$genestrand $LCY$readstrand $convtype $LPR$flag$N\t$LGN$data{$nopkFile}{total}$N\t$nopkFile\n\n";
}


my @strand = qw(Pos Neg);
my @convtype = qw(CG CH GC GH);
my @strandconvtype = qw(Pos_CG Neg_GC Pos_GC Neg_CG Pos_CH Neg_GH Pos_GH Neg_CH);
open (my $out1, ">", "$footPeakFolder/99_FOOTSTATS/1_STATS.txt") or die;
print $out1 "BC\tplasmid\tdesc\ttx";
#foreach my $strand (@strandtrand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
foreach my $strand (@strand[0..@strand-1]) {
	print $out1 "\t$strand";
}
#	}
#}
#foreach my $strand (@strand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\ttotal.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tpeakmean.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tnopkmean.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tpeak25.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tnopk25.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tpeak50.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tnopk50.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tpeak75.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tnopk75.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tpeak100.$strandconvtype";
}
foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
	print $out1 "\tnopk100.$strandconvtype";
}
#	}
#}
print $out1 "\n";
foreach my $bc (sort keys %final) {
	foreach my $plasmiddesc (sort keys %{$final{$bc}}) {
		my ($plasmid, $desc) = $plasmiddesc =~ /^PLASMID(.+);DESC(.+)$/;
		my ($tx) = $desc =~ /(TX_RH1|TX|NT|NT_RH1)$/;
		my ($desc2) = $desc =~ /^(.+)_(TX_RH1|TX|NT|NT_RH1)$/;
		$tx = "NA" if not defined $tx;
		$desc2 = $desc if not defined $desc2;
		print $out1 "$bc\t$plasmid\t$desc2\t$tx";
		my $peakpercprint = "";
		my $totalreadprint = "";
		my $peakmean = "";
		my $nopkmean = "";
		my $peak25print = "";
		my $peak50print = "";
		my $peak75print = "";
		my $peak100print = "";
		my $nopk25print = "";
		my $nopk50print = "";
		my $nopk75print = "";
		my $nopk100print = "";
		foreach my $strand (@strand[0..@strand-1]) {	
			my $convtype = "CH";
			my $totalread = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread};
			$totalread = 0 if not defined $totalread;
			$totalreadprint .= "\t$totalread";
		}
		foreach my $strandconvtype (@strandconvtype[0..@strandconvtype-1]) {
		#foreach my $strand (@strand[0..@strand-1]) {#sort keys %{$final{$bc}{$plasmiddesc}}) {
		#	foreach my $convtype (@convtype[0..@convtype-1]) {#sort keys %{$final{$bc}{$plasmiddesc}{$strand}}) {
				my ($strand, $convtype) = $strandconvtype =~ /^(\w+)_(\w+)$/;
				my $flag = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{flag};				
				my $peak = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak};				
				my $nopk = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk};				
				my $totalread = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{totalread};				
				my $peakperc = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakperc};				
				my $nopkperc = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkperc};
				my $peakmean1 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peakmean}; $peakmean1 = 0 if not defined $peakmean1;
				my $nopkmean1 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopkmean}; $nopkmean1 = 0 if not defined $nopkmean1;
				my $peakprint25 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak25}; $peakprint25 = 0 if not defined $peakprint25;
				my $peakprint50 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak50}; $peakprint50 = 0 if not defined $peakprint50;
				my $peakprint75 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak75}; $peakprint75 = 0 if not defined $peakprint75;
				my $peakprint100 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{peak100}; $peakprint100 = 0 if not defined $peakprint100;
				my $nopkprint25 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk25}; $nopkprint25 = 0 if not defined $nopkprint25;
				my $nopkprint50 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk50}; $nopkprint50 = 0 if not defined $nopkprint50;
				my $nopkprint75 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk75}; $nopkprint75 = 0 if not defined $nopkprint75;
				my $nopkprint100 = $final{$bc}{$plasmiddesc}{$strand}{$convtype}{nopk100}; $nopkprint100 = 0 if not defined $nopkprint100;
				$totalread = 0 if not defined $totalread;
				$peakperc = 0 if not defined $peakperc;
				$nopkperc = 0 if not defined $nopkperc;
				$peakmean .= "\t$peakmean1";
				$nopkmean .= "\t$nopkmean1";
				$peak25print .=  "\t$peakprint25";
				$peak50print .=  "\t$peakprint50";
				$peak75print .=  "\t$peakprint75";
				$peak100print .=  "\t$peakprint100";
				$nopk25print .=  "\t$nopkprint25";
				$nopk50print .=  "\t$nopkprint50";
				$nopk75print .=  "\t$nopkprint75";
				$nopk100print .=  "\t$nopkprint100";
				$peakpercprint .= "\t$peakperc";
#				$nopkprint .=  "\t$nopkprint1";
		#	}
		#}
		#	print "$bc\t$plasmid\t$desc2\t$strandconvtype\t$peakmean1\t$nopkmean1\t$peakprint25\t$peakprint50\t$peakprint100\n\n";
		}

		print $out1 "$totalreadprint";
		print $out1 "$peakpercprint";
		print $out1 "$peakmean";
		print $out1 "$nopkmean";
		print $out1 "$peak25print";
		print $out1 "$nopk25print";
		print $out1 "$peak50print";
		print $out1 "$nopk50print";
		print $out1 "$peak75print";
		print $out1 "$nopk75print";
		print $out1 "$peak100print";
		print $out1 "$nopk100print";
		#print $out1 "$nopkprint";
		print $out1 "\n";
	}
}

print "SUCCESS!!!\n";
exit 0;

sub check0 {
	my ($df, $ops) = @_;
	my @ops = @{$ops};
	my @ops0 = @{$ops[0]};
	my ($temphash2, $totaltemp2);
	if (not defined $df) {
	#	print "[" . join(" ", @{$ops[$0]}) . "]: ${LGN}0$N\n";
		return($df, 0);
	}
	else {
	#	print join(" ", @ops0) . "\n\n";
	}
	($temphash2, $totaltemp2) = check($df, @ops0);
	#print "[" . join(" ", @ops0) . "]: $LGN$totaltemp2$N\n";
	for (my $i = 1; $i < @ops; $i++) {
		($temphash2, $totaltemp2) = check($temphash2, @{$ops[$i]});
		#print "[" . join(" ", @{$ops[$i]}) . "]: $LGN$totaltemp2$N\n";
	}
	return($temphash2, $totaltemp2);
}
sub check {
	my ($df, $colname, $op, $num2) = @_;
	#print "$colname $op $num2\n";
	if (not defined $df) {
		return($df, 0);
	}
	if (keys (%{$df}) == 0) {
		return($df, 0);
	}
	#my @dfarr = @{$dfarr{$colname}};
	my %df = %{$df};
	my %df2;
	my %total;
	for (my $i = 0; $i < @{$df{i}}; $i++) {
		my $val = $df{$colname}[$i];
		next if not defined $val;
		next if $val !~ /^\-?\d+\.?\d*e?\-?\d*$/;
		if ($op eq "gt") {
			next if $val <= $num2;
		}
		elsif ($op eq "gte") {
			next if $val < $num2;
		}
		elsif ($op eq "lt") {
			next if $val >= $num2;
		}
		elsif ($op eq "lte") {
			next if $val > $num2;
		}
		elsif ($op eq "eq") {
			next if $val != $num2;
		}
		my $readname = $df{read}[$i];
		$total{$readname} ++;
		foreach my $key (keys %df) {
			print "$key $df{$key}[$i]\n" if not defined $readname;
			push(@{$df2{$key}}, $df{$key}[$i]);
		}
		die "Undef read\n" if not defined $readname;
	}
	my $total2 = (keys %total);
	return(\%df2, $total2);
}
sub footStats_parse_footPeak_logFile {
	my ($footPeakFolder, $footPeak_logFile, $outLog) = @_;
	my %opts; my $debugprint = "\n\n$YW<--------- 0. PARSING LOGFILES ----------$N\n\n"; my $optprint = "";
	open (my $footPeak_logFileIn, "<", $footPeak_logFile) or DIELOG($outLog, "Cannot read from $footPeak_logFile: $!\n");
	while (my $line = <$footPeak_logFileIn>) {
		chomp($line);
		my $check = 0;
		if ($line =~ />Run Params$/) {
				$check ++;
				die "0: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 1;
				while ($line !~ />Options:$/) {
						my ($param, $value) = $line =~ /^([\w ]+[a-zA-Z]+)[ \t]+:[ \t]+(.+)$/;
						$line =~ />Run Params$/ ? $debugprint .= "\n$LCY 0.$N footPeak $line:\n" : (defined $param and $param !~ /^Run script/) ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
						$opts{footPeak}{$param} = $value if (defined $param and defined $value and $param !~ /^Run script/);
						if (defined $param and $param eq "Run script short") {
								my @values = split(" -", $value); shift(@values);
								foreach my $values (@values) {
										my ($param2, $value2) = $values =~ /^(\w) (.+)$/;
										$debugprint .= "footLoopFolder = $LCY$value2$N\n" if $values =~ /^n /;
										$opts{footLoop2}{n} = $value2 if $values =~ /^n /;
										$opts{footPeak}{$param2} = $value2;
								}
						}
						$line = <$footPeak_logFileIn>; chomp($line);
				}
				$optprint .= ">footPeak\n";
				foreach my $param (sort keys %{$opts{footPeak}}) {
						$optprint .= "$param=$opts{footPeak}{$param}\n";
				}
		}
		if ($line =~ /^>Options/) {
				$check ++;
				die "1: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 2;
				while ($line !~ />Run Params from footLoop.pl/) {
						my ($param, $desc, $value) = $line =~ /^\-(\w)\s*(\w+)\s*:\s*([a-zA-Z0-9]+)$/;
						   ($param, $value) = $line =~ /^\-(\w)\s*:\s*([a-zA-Z0-9]+)$/ if not defined $param;
						$line =~ />Options/ ? $debugprint .= "\n$LGN 1.$N footPeak $line: " : defined $param ? $debugprint .= "$LCY$param$N=$value;" : $debugprint .= "";
						$opts{footPeak2}{$param} = $value if (defined $param and defined $value);
						$line = <$footPeak_logFileIn>; chomp($line);
				}
				foreach my $param (sort keys %{$opts{footPeak2}}) {
						$optprint .= "$param=$opts{footPeak2}{$param}\n";
				}
		}
		if ($line =~ />Run Params from footLoop.pl/) {
				$check ++; my $debugprint2;
				die "2: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 3;
				while ($line !~ /^>Options from footLoop.pl/) {
						my ($param, $value) = $line =~ /^footLoop\s*(\w+)\s*:[ \t]+(.+)$/;
						if ($line =~ /from footLoop.pl logfile/) {($param, $value) = $line =~ /(logfile)=\.?\/?(.+)$/; $value = $opts{footLoop2}{n} . $value; $debugprint2 .= "\t- $LCY$param$N=$value\n";}
						$line =~ />Run Params/ ? $debugprint .= "\n$LCY 2.$N footLoop $line:\n$debugprint2" : defined $param ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
						$opts{footLoop}{$param} = $value if (defined $param and defined $value and $param ne "origDir");
						$line = <$footPeak_logFileIn>; chomp($line);
				}
				$optprint .= ">footLoop\n";
				foreach my $param (sort keys %{$opts{footLoop}}) {
						$optprint .= "$param=$opts{footLoop}{$param}\n";
				}

		}
		if ($line =~ /^>Options from footLoop.pl/) {
				$check ++;
				die "3: ${LRD}ERROR!$N footPeak.pl logfile $LCY$footPeak_logFile$N is corrupted or different version than 2.95!\n" unless $check == 4;
				while ($line !~ /^.+[\-]+\>/) {
						my ($param, $value) = $line =~ /^\-(\w)\s*:\s+(.+)$/;
						$opts{footLoop2}{$param} = $value if (defined $param and defined $value);
						$line =~ />Options from footLoop/ ? $debugprint .= "\n$YW 3.$N footLoop $line:\n" : defined $param ? $debugprint .= "\t- $LCY$param$N=$value\n" : $debugprint .= "";
						$line = <$footPeak_logFileIn>; chomp($line);
				}
				foreach my $param (sort keys %{$opts{footLoop2}}) {
						$optprint .= "$param=$opts{footLoop2}{$param}\n";
				}
		}
	#	   my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $bedLine);
	}
	$debugprint .= "\n\n$YW------------------------------------->$N\n\n\n";
	makedir("$footPeakFolder/99_FOOTSTATS");
	open (my $out, ">", "$footPeakFolder/99_FOOTSTATS/.PARAMS") or DIELOG($outLog, "Failed to write to $footPeakFolder/99_FOOTSTATS/.PARAMS: $!\n");
	print $out "$optprint\n";
	close $out;
	LOG($outLog, $debugprint,"NA");
	return \%opts;
}
sub parse_geneIndexFile {
   my ($geneIndexFile, $outLog) = @_;
   my %coor;
   LOG($outLog, "${LCY}geneIndexFile$N=$geneIndexFile\n","NA");
   die "geneindexFile does not exist!\n" if not defined $geneIndexFile;
   open (my $in, "<", $geneIndexFile) or DIELOG($outLog, "Failed to read from $geneIndexFile: $!\n");
   while (my $line = <$in>) {
	  chomp($line);
	  my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	  $gene = uc($gene);
	  $coor{uc($gene)}{chr} = $chr;
	  $coor{uc($gene)}{beg} = $beg;
	  $coor{uc($gene)}{end} = $end;
	  $coor{uc($gene)}{strand} = $strand;
	  #print "chr=$chr:$beg-$end gene=$gene, strand=$strand\n";
   }
   close $in;
   return \%coor;
}

sub check_sanity {

	my $usage = "

-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N [Optional: -G <genewant>] -n $LCY<footPeak_folder>$N

";

	die $usage unless defined $opt_n and -d $opt_n;

	my ($footPeakFolder) = $opt_n;
	my ($logFile) = $footPeakFolder . "/footStats_logFile.txt";
	my ($outDir) = $footPeakFolder . "/99_FOOTSTATS/";
	my $optG = defined $opt_G ? "-G $opt_G " : "";
	makedir("$footPeakFolder/99_FOOTSTATS/") if not -d "$footPeakFolder/99_FOOTSTATS/";
	open (my $outLog, ">", $logFile) or die "\nFailed to write to $LCY$logFile$N: $!\n";
	print $outLog "$0 version $version\n\nRun script: $0 $optG-n $opt_n\n";
	return($footPeakFolder, $outLog);
}



=comment
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Processing BAMFile$N $LCY$BAMFile$N\n");
my (%BAMData) = %{parse_BAMFile($BAMFile, $outLog)};
foreach my $gene (sort keys %BAMData) {
	my @key = qw(Pos Neg);
	my $total = $BAMData{uc($gene)}{total};
	my $pos = $BAMData{uc($gene)}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $BAMData{uc($gene)}{Neg}; $neg = 0 if not defined $neg;
#	print "$gene\t$total\t$pos\t$neg\n";
}

my ($BAMFixedFile) = $opts->{footLoop}{BAMFixed};
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Processing BamFixedFile$N $LCY$BAMFixedFile$N\n");
my ($origDir) = getFilename($BAMFixedFile, "folder");
my (%BAMFixedData) = %{parse_BAMFile($BAMFixedFile, $outLog)};

LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "${LPR}Printing each gene in BamFixedData$N\n");
my %dataorig;
my %print0;
$print0{header} = "gene\tsam_uread\tsam_uread_pos\tsam_uread_neg\tsamfix_uread_pos\tsamfix_uread_neg\tsamfix_pos_pos\tsamfix_neg_neg\tsamfix_pos_neg\tsamfix_neg_pos\tsamfinal_uread_pos\tsamfinal_uread_neg\tsamfinal_uread_unk\tinclude_cpg";
foreach my $gene (sort keys %BAMFixedData) {
	my @key = qw(Pos Neg);
	my $total = $BAMFixedData{uc($gene)}{total};
	my $sampos = $BAMData{uc($gene)}{Pos}; $sampos = 0 if not defined $sampos;
	my $samneg = $BAMData{uc($gene)}{Neg}; $samneg = 0 if not defined $samneg;
	my $pos = $BAMFixedData{uc($gene)}{Pos}; $pos = 0 if not defined $pos;
	my $neg = $BAMFixedData{uc($gene)}{Neg}; $neg = 0 if not defined $neg;
	die "Discrepancy betweensam and fixed sam pos and neg in gene=$gene (pos=$pos, sampos=$sampos)\n" if $sampos ne $pos;
	die "Discrepancy betweensam and fixed sam neg and neg in gene=$gene (neg=$neg, samneg=$samneg)\n" if $samneg ne $neg;
	my $pos2 = $BAMFixedData{uc($gene)}{strand}{Pos}; $pos2 = 0 if not defined $pos2;
	my $neg2 = $BAMFixedData{uc($gene)}{strand}{Neg}; $neg2 = 0 if not defined $neg2;
	my $posneg = $BAMFixedData{uc($gene)}{change}{Pos}{Neg}; $posneg = 0 if not defined $posneg;
	my $negpos = $BAMFixedData{uc($gene)}{change}{Neg}{Pos}; $negpos = 0 if not defined $negpos;
	my $pospos = $BAMFixedData{uc($gene)}{change}{Pos}{Pos}; $pospos = 0 if not defined $pospos;
	my $negneg = $BAMFixedData{uc($gene)}{change}{Neg}{Neg}; $negneg = 0 if not defined $negneg;
	my $posorig = "$origDir/$gene\_Pos.orig"; ($posorig) = -e $posorig ? `wc -l $posorig` =~ /^\s*(\d+)/ : 0;
	my $negorig = "$origDir/$gene\_Neg.orig"; ($negorig) = -e $negorig ? `wc -l $negorig` =~ /^\s*(\d+)/ : 0;
	my $unkorig = "$origDir/$gene\_Unk.orig"; ($unkorig) = -e $unkorig ? `wc -l $unkorig` =~ /^\s*(\d+)/ : 0;
	$dataorig{uc($gene)}{posorig} = $posorig;
	$dataorig{uc($gene)}{negorig} = $negorig;
	$dataorig{uc($gene)}{unkorig} = $unkorig;
	my $include_cpg = $BAMData{uc($gene)}{include_cpg};
	my @print0curr = ($gene,$total,$pos,$neg,$pos2,$neg2,$pospos,$negneg,$posneg,$negpos,$posorig,$negorig,$unkorig,$include_cpg);
	for (my $i = 0; $i < @print0curr; $i++) {
		$print0curr[$i] = "NA" if not defined $print0curr[$i]; 
		$print0curr[$i] = "NA" if $print0curr[$i] =~ /^\s*$/;
	}
	$print0{gene}{uc($gene)} = join("\t", @print0curr);
	
}


sub parse_BAMFile {
	my ($BAMFile, $outLog) = @_;
	my %data;
	my $in;
	if ($BAMFile =~ /.gz$/) {
		open ($in, "zcat < $BAMFile|") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	elsif ($BAMFile =~ /.sam$/) {
		open ($in, "<", "$BAMFile") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	else {
		open ($in, "samtools view $BAMFile|") or DIELOG($outLog, date() . " Failed to read from BAMFile $LCY$BAMFile$N: $!\n");
	}
	my $include_cpg = $BAMFile =~ /PCB\d+_bcBC\d+_plasmid/ ? "TRUE" : "FALSE";
	while (my $line = <$in>) {
		chomp($line);
		my @arr = split("\t", $line);
		next if @arr < 6;
		my ($read, $flag, $gene, $pos, $qual) = @arr;
		my ($flag2, $type);
		if ($arr[1] !~ /^(0|16)$/) {
			($read, $type, $flag, $flag2, $gene) = @arr;
		}
		my $strand = $flag eq 16 ? "Neg" : $flag eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$BAMFile$N, flag=$flag, line:\n\n$line\n\n");
		$data{uc($gene)}{$strand} ++;
		$data{uc($gene)}{total} ++;
		$data{uc($gene)}{include_cpg} = $include_cpg;
		if ($arr[1] !~ /^(0|16)$/) {
			my $strand2 = $flag2 eq 16 ? "Neg" : $flag2 eq 0 ? "Pos" : DIELOG($outLog, "Failed to parse strand from sam file=$LCY$BAMFile$N, flag2=$flag2, line:\n\n$line\n\n");
			$data{uc($gene)}{strand}{$strand2} ++;
			$data{uc($gene)}{change}{$strand}{$strand2} ++;
		}
	}
	foreach my $gene (sort keys %data) {
		foreach my $strand2 (sort keys %{$data{uc($gene)}{strand}}) {
			LOG($outLog, date() . "${LPR}gene info$N: gene=$gene strand=$strand2 total=$data{uc($gene)}{strand}{$strand2}\n","NA");
		}
	}
	close $in;
	return \%data;
}

# get fixed sam file
my ($fixedBAMFile) = $opts->{footLoop}{BAMFixed};



my (@TXTFile) = <$footPeakFolder/99_FOOTSTATS/.PEAKSTATSTEMP/.0_*.TXT>;
if (@TXTFile == 0) {@TXTFile = <$footPeakFolder/.0_*.TXT>;}

my $data;
my $TXTFileInd = 0;
my $totalTXTFile = @TXTFile;
foreach my $TXTFile (@TXTFile[0..@TXTFile-1]) {
	my ($folder1, $fileName1) = getFilename($TXTFile, "folderfull");
	my ($footName) = $fileName1 =~ /.0_RESULTS_(.+).TXT$/;
	if ($TXTFileInd <= 20 or ($TXTFileInd > 20 and $TXTFileInd % 100 == 0)) {
		LOG($outLog, date() . "${LPR}Parsing TXTFile$N $LGN$TXTFileInd/$totalTXTFile$N $LCY$fileName1$N footName = $LCY$footName$N\n");
	}
	else {
		LOG($outLog, date() . "${LPR}Parsing TXTFile$N $LGN$TXTFileInd/$totalTXTFile$N $LCY$fileName1$N footName = $LCY$footName$N\n","NA");
	}
	my $parseName = parseName($footName . "_CH");
	my @header = qw(label gene strand window thres conv pcb bc plasmid desc);
	my ($label, $gene, $strand, $window, $thres, $conv, $pcb, $bc, $plasmid, $desc) = @{$parseName->{array}};
	$conv = "";
	$data = parseTXT($data, $TXTFile, $label, $gene, $strand, $window, $thres);
	$TXTFileInd ++;
}

my %print1;
my %print;
my $LABEL;

my $outPeakStatsFile = "$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT";
LOG($outLog, date() . "\n------------------------\n");
LOG($outLog, date() . "\n${LPR}Printing each gene into $LCY$outPeakStatsFile$N\n");
open (my $out, ">", $outPeakStatsFile) or DIELOG($outLog, date() . " Failed to write to $outPeakStatsFile: $!\n");
print $out "type\tlabel\tgene\tread_strand\twindow\tthres\tconv\tread_unique_total\tread_unique_peak_total\tread_unique_peak_fraction\tfootPeakFolder\tpeakFile\tinclude_cpg\n";
my $outPeakStatsInd = 0;
foreach my $label (sort keys %{$data}) {
	foreach my $gene (sort keys %{$data->{$label}}) {
		foreach my $read_strand (sort keys %{$data->{$label}{uc($gene)}}) {
			foreach my $window (sort keys %{$data->{$label}{uc($gene)}{$read_strand}}) {
				foreach my $thres (sort keys %{$data->{$label}{uc($gene)}{$read_strand}{$window}}) {
					foreach my $conv (sort keys %{$data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}}) {
						my $strand = $coor->{uc($gene)}{strand}; $strand = "Pos" if $strand eq "+"; $strand = "Neg" if $strand eq "-";
						my $goodconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "CG" : $read_strand eq "Neg" ? "GC" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						$goodconv = $goodconv ne "H" ? $goodconv : $read_strand eq "Pos" ? "CH" : $read_strand eq "Neg" ? "GH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get goodconv on label=$label strand eq $read_strand\n");
						my $goodconv2 = $goodconv eq "CH" ? "CG" : $goodconv eq "GH" ? "GC" : $goodconv eq "CG" ? "CH" : "GH";
						my $revconv = $label =~ /^PCB\d+$/ ? "H" : $read_strand eq "Pos" ? "GC" : $read_strand eq "Neg" ? "CG" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						$revconv = $revconv ne "H" ? $revconv : $read_strand eq "Pos" ? "GH" : $read_strand eq "Neg" ? "CH" : $read_strand eq "Unk" ? "UNK" : DIELOG($outLog, "Failed to get revconv on label=$label strand eq $read_strand\n");
						my $type;
						if ($outPeakStatsInd <= 100 or ($outPeakStatsInd > 100 and $outPeakStatsInd % 1000 == 0)) {
							LOG($outLog, date() . "Done $LGN$outPeakStatsInd$N: label=$label gene=$gene read_strand = $read_strand gene_strand = $strand goodconv = $goodconv\n");
						}
						else {
							LOG($outLog, date() . "Done $LGN$outPeakStatsInd$N: label=$label gene=$gene read_strand = $read_strand gene_strand = $strand goodconv = $goodconv\n","NA");
						}
						$outPeakStatsInd ++;

						my $include_cpg = $label =~ /^PCB\d+$/ ? "FALSE" : "TRUE";
						#my $cpg2 = $label =~ /^PCB\d+$/ ? "_C" : "";
						$type = 	$read_strand eq "UNK" ? "UNK" :
									($strand eq $read_strand and $conv eq $goodconv) ? "PEAK" :
									($strand ne $read_strand and $conv eq $goodconv) ? "PEAK_TEMP" : 
									($strand eq $read_strand and $conv eq $revconv) ? "PEAK_RCONV" :
									($strand ne $read_strand and $conv eq $revconv) ? "PEAK_TEMP_RCONV" :
									($strand eq $read_strand and $conv eq $goodconv2) ? "PEAK_OTHERC" :
									($strand ne $read_strand and $conv eq $goodconv2) ? "PEAK_TEMP_OTHERC" : "OTHERS";
						my ($currfootPeakFolder) = getFullpath($data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{footPeakFolder});
						my ($currpeakFile) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{peakFile};
						my ($read_unique_total) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_total};
						my ($read_unique_peak_total) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_total};
						my ($read_unique_peak_fraction) = $data->{$label}{uc($gene)}{$read_strand}{$window}{$thres}{$conv}{read_unique_peak_fraction} / 100;
						$print1{uc($gene)}{$type} += $read_unique_peak_total;
						$print1{uc($gene)}{strand} = $strand;
						$print1{uc($gene)}{OTHERS_TOTAL} ++ if $type eq "OTHERS";
						
						$print{$type} .= "$type\t$label\t$gene\t$read_strand\t$window\t$thres\t$conv\t$read_unique_total\t$read_unique_peak_total\t$read_unique_peak_fraction\t$currfootPeakFolder\t$currpeakFile\t$include_cpg\n";
						$print{ALL} = "$type\t$label\t$gene\t$read_strand\t$window\t$thres\t$conv\t0\t0\t0\t$currfootPeakFolder\t$currpeakFile\t$include_cpg\n";
						$LABEL = $label;
					}
				}
			}
		}
	}
}

open (my $out2, ">", "$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT") or die;
print $out2 "label\tstrand\t$print0{header}";
my @types = qw(PEAK PEAK_TEMP PEAK_RCONV PEAK_TEMP_RCONV OTHERS); 
for (my $i = 0; $i < @types; $i++) { 
	if (not defined $print{$types[$i]}) {
		$print{$types[$i]} = $print{ALL};
	}
	print $out $print{$types[$i]};
	print $out2 "\t$types[$i]\t$types[$i]/total";
}
print $out2 "\n";

foreach my $gene (sort keys %{$print0{gene}}) {
	my $print0 = $print0{gene}{uc($gene)};
	my $strand2 = $coor->{uc($gene)}{strand}; 
	my $strand = $strand2 eq "+" ? "Pos" : $strand2 eq "-" ? "Neg" : "UNKSTRAND";
	die "Undefined strand (strand=$strand, strand2=$strand2) at gne=$gene\n" if $strand eq "UNKSTRAND";
	my $posorig = $dataorig{uc($gene)}{posorig};
	my $negorig = $dataorig{uc($gene)}{negorig};
	my $unkorig = $dataorig{uc($gene)}{unkorig};
	print $out2 "$LABEL\t$strand\t$print0";
	for (my $i = 0; $i < @types-1; $i++) { 
		my $print1 = $print1{uc($gene)}{$types[$i]}; $print1 = 0 if not defined $print1;
		my $print1denomPos = $strand eq "Pos" ? $posorig : $strand eq "Neg" ? $negorig : $unkorig; $print1denomPos = 1 if $print1denomPos == 0;
		my $print1denomNeg = $strand eq "Pos" ? $negorig : $strand eq "Neg" ? $posorig : $unkorig; $print1denomNeg = 1 if $print1denomNeg == 0;
	#	print "type=$types[$i] strand=$strand print1=$print1 denom = $print1denomPos neg = $print1denomNeg\n";
		print $out2 "\t$print1";
		if ($i < @types - 1) {
			print $out2 "\t" . int(100 * $print1 / $print1denomPos +0.5)/100 if $i % 2 == 0;
			print $out2 "\t" . int(100 * $print1 / $print1denomNeg +0.5)/100 if $i % 2 == 1;
			LOG($outLog, "strand=$strand, $print1denomPos\n","NA") if $i % 2 == 0;
			LOG($outLog, "strand=$strand, $print1denomNeg\n","NA") if $i % 2 == 1;
		}
	}
	my $printOTHERS = $print1{uc($gene)}{OTHERS}; $printOTHERS = 0 if not defined $printOTHERS;
	my $printOTHERStot = $print1{uc($gene)}{OTHERS_TOTAL}; $printOTHERStot = 1 if not defined $printOTHERStot;
	
	my $printOTHERSperc = int($printOTHERS / $printOTHERStot * 10)/10;
	print $out2 "\t$printOTHERS\t$printOTHERSperc\n";
}
LOG($outLog, date() . "$LCY$opt_n$N Done!\n\nOutput:
$footPeakFolder/99_FOOTSTATS/1_PEAKSTATS.TXT
$footPeakFolder/99_FOOTSTATS/0_SUMMARY.TXT\n");


sub parseTXT {
	my ($data, $TXTFile, $label, $gene, $strand, $window, $thres) = @_;
	my @header = qw(footPeakFolder peakFile gene conv read_unique_total read_unique_peak_total read_unique_peak_fraction peakType);

	open (my $in1, "<", $TXTFile) or die "Cannot read from $TXTFile: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		if ($line =~ /^#/) {
			next;
		}
		my @arr = split("\t", $line);
		if (@header != @arr) {
			LOG($outLog, "total number of header isn't the same as total column in file $LCY$TXTFile$N!\n\n$LPR$line$N\n\n");
			my $max = @header > @arr ? @header : @arr;
			for (my $i = 0; $i < $max; $i++) {
				my $headerVal = defined $header[$i] ? $header[$i] : "UNDEF";
				my $arrVal = defined $arr[$i] ? $arr[$i] : "UNDEF";
				print "$i header=$headerVal val=$arrVal\n";
			}
			DIELOG($outLog, "\n");
		}
#		print "\n$line\n";
		my $conv = $arr[3];
		for (my $i = 0; $i < @arr; $i++) {
			next if $i == 3 or $i == 2;
			$data->{$label}{uc($gene)}{$strand}{$window}{$thres}{$conv}{$header[$i]} = $arr[$i];
		}
	}
	close $in1;
	return($data);
}

=cut
__END__
sub parse_geneIndexFile {
   my ($footLoopFolder, $outLog) = @_;
   my ($geneIndexFile) = <$footLoopFolder/*.bed>;
   my %coor;
   LOG($outLog, "${LCY}geneIndexFile$N=$geneIndexFile\n","NA");
   die "geneindexFile does not exist!\n" if not defined $geneIndexFile;
   open (my $in, "<", $geneIndexFile) or DIELOG($outLog, "Failed to read from $geneIndexFile: $!\n");
   while (my $line = <$in>) {
	  chomp($line);
	  my ($chr, $beg, $end, $gene, $zero, $strand) = split("\t", $line);
	  $gene = uc($gene);
	  $coor{uc($gene)}{chr} = $chr;
	  $coor{uc($gene)}{beg} = $beg;
	  $coor{uc($gene)}{end} = $end;
	  $coor{uc($gene)}{strand} = $strand;
	  #print "chr=$chr:$beg-$end gene=$gene, strand=$strand\n";
   }
   close $in;
   return \%coor;
}

__END__
my %data;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

170802_pcb05_RnaseHneg_ccs_3minFP.fastq.gz.rmdup.fq.gz_PEAK_20W35T	  CALM3_Neg_20_0.35_GC.PEAK	   CALM3   GC	  160	 3	   1.9


