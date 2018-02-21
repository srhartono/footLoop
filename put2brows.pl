#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_b $opt_f $opt_n $opt_l $opt_d);
getopts("vb:f:n:l:d:");

BEGIN {
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";


my ($footFolder, $peakFolder, $dupeFolder, $bedFile) = ($opt_f, $opt_n, $opt_d, $opt_b);
die "\nusage: $YW$0$N -f $CY<footPeak -n output dir>$N -d <dupeFolder> $LGN [optional: -b bedFile index or -f footLoop folder]$N\n\n" unless defined $opt_n and -d $opt_n;
die "\nPlease supply -d dupefolder\n" if not defined $opt_d or not -d $opt_d;
die "\nPlease either supply -b (bedFile index in footloop folder) or -f (footloop folder)>$N\n\n" if not defined $opt_b and not defined $opt_f;

$bedFile = $opt_b if defined $opt_b;
my %coor;
if (not defined $opt_b) {
	my (@bedFile) = <$footFolder/*.bed>;
	die "multiple bedFiles in $footFolder please use -b with one of them!\n" if @bedFile > 1;
	die "no bedFile index in $footFolder! Please manually supply one with -b\n" if @bedFile == 0;
	$bedFile = $bedFile[0];
}
my @line = `cat $bedFile`;
foreach my $line (@line) {
	chomp($line);
	my ($chr, $beg, $end, $gene) = split("\t", $line);
	$gene = uc($gene);
	$beg = 1 if $beg == 0;
	$coor{$gene}{end} = $end;
	$coor{$gene}{beg} = $beg;
	$coor{$gene}{chr} = $chr;
	print "$gene: $beg\n";
}

my $label = defined $opt_l ? $opt_l : $peakFolder;
$label =~ s/\//_/g; $label =~ s/_$//g;
print "label=$label\n";

# dupe
open (my $out1, ">", "$peakFolder/allgtf.gtf") or die "Cannot write to $peakFolder/allgtf.gtf; $!\n";
foreach my $genez (sort keys %coor) {
print "Doing $genez\n";
my %dupe; my $maxgrp = 0;
my @dupe = <$dupeFolder/*$genez*.order>;
if (@dupe == 0) {print "$genez: There is no .order file in $LCY$dupeFolder$N\n" and next;}
elsif (@dupe > 1) {print "$genez: There are more than 1 .order files in $LCY$dupeFolder$N\n" and next;}
else {
	@dupe = `cut -f1-3 $dupe[0]`;
	foreach my $dupe (@dupe) {
		chomp($dupe);
		my ($read, $grp, $id) = split("\t", $dupe);
		push(@{$dupe{grp}{$grp}}, $read);
		$dupe{read}{$read} = $grp;
		$maxgrp = $grp if $maxgrp < $grp;
	}
}
my @calls = <$peakFolder/.CALL/*$genez*.PEAK.out>;
print "- $LGN" . join("$N\n- $LGN", @calls) . "$N\n\n";
foreach my $input1 (@calls) {
	my $fullname = getFullpath($input1);
	my ($name1) = $fullname =~ /.+(PCB\d+)/; $name1 = "" if not defined $name1;
	my ($name2) = $fullname =~ /.+_ccs_(\w+).fastq./ if $fullname =~ /_ccs_BC/; $name2 = "" if not defined $name2;
	$label = "$name1\_$name2" if $name1 ne "" and $name2 ne "" and not defined $opt_l;
	my ($fileName1) = getFilename($input1,'full');
	$fileName1 =~ s/.PEAK.out//;
	my ($gene) = $fileName1 =~ /^(.+)_(Pos|Neg|Unk)/;
	die "gene undefined filename = $fileName1, input1=$input1\n" if not defined $gene;
	$gene = uc($gene);
	my $end = $coor{$gene}{end};
	my $beg = $coor{$gene}{beg};
	my $chr = $coor{$gene}{chr};
	die "gene $gene cannot find beg in bed file $bedFile\n" if not defined $beg;
	my %data;
	my %temp;
	print "Doing $input1\n";
	my ($strand) = $input1 =~ /Pos/ ? "+" : $input1 =~ /Neg/ ? "-" : $input1 =~ /Unk/ ? "+" : die "Failed to get strand from file=$input1: $!\n";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	my $color = $input1 =~ /Pos/ ? "200,50,0" : $input1 =~ /Neg/ ? "0,50,200" : "0,155,0";
	my %read;
	my $track = "track type=gtf name=$label\_$fileName1 color=$color visibility=0\n";
	my $lastgrp = -1; my $group2 = 0;
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		my ($read, @val) = split("\t", $line);
		my $group = $dupe{read}{$read}; die "Dead at gene=$gene read=$read not defined group!\n" if not defined $group;
		my $vals = join("", @val);
		# 001100
		# 012345
		# 2 and 4
		# i = 2; i < 4; i++
		my ($seqborder0) = $vals =~ /^(0*)[1-9]/; 
		$seqborder0 = defined $seqborder0 ? length($seqborder0) + $beg : 0 + $beg;
		my ($seqborder1) = $vals =~ /^(.+[1-9])0*$/; 
		$seqborder1 = defined $seqborder1 ? length($seqborder1) + $beg : @val + $beg;
#		my $group2 = int($group / $maxgrp * 500+0.5) + 1;
		my $begGTF;
		if ($lastgrp eq -1 or ($lastgrp ne -1 and $lastgrp ne $group)) {
#			$group2 ++;
      	$begGTF  = "$chr\tNA\tgene\tMUHGROUP2\t$seqborder1\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2\"\n";
         $begGTF .= "$chr\tNA\texon\tMUHGROUP2\t$end\t0\t$strand\t0\tgene_id \"MUHGROUP2.2\"; transcript_id \"MUHGROUP2.2\"\n";
		}
		$lastgrp = $group;
			$begGTF .= "$chr\tNA\ttranscript\tMUHGROUP2\t$seqborder1\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
         $begGTF .= "$chr\tNA\texon\tMUHGROUP2\tMUHGROUP2\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
  #       $begGTF .= "$chr\tNA\texon\t$seqborder0\t$seqborder0\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"";
#         $begGTF .= "$chr\tNA\texon\t$seqborder0\t$seqborder0\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"";
      my $endGTF  = "$chr\tNA\texon\t$seqborder1\t$seqborder1\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
         $endGTF .= "$chr\tNA\texon\t$end\t$end\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
#         $endGTF .= "$chr\tNA\texon\t$seqborder1\t$seqborder1\t0\t$strand\t0\tgene_id \"MUHGROUP2.2\"; transcript_id \"MUHGROUP2.-2\"";
		my $midGTF  = "\n";
		my ($currbeg, $currend) = (-1,-1);
		for (my $i = $seqborder0; $i < $seqborder1; $i++) {
			my $val = $val[$i];
			my $pos = $i + $beg;
			next if $val !~ /[89]/;
			if ($currbeg == -1) {
				$currbeg = $pos;
			}
			elsif ($pos != $currend + 1) {
				$midGTF .= "$chr\tNA\texon\t$currbeg\t$currend\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
				$currbeg = $pos; 
			}
			$currend = $pos;
			if ($i == $seqborder1 - 1) {
				$midGTF .= "$chr\tNA\texon\t$currbeg\t$currend\t0\t$strand\t0\tgene_id \"MUHGROUP2\"; transcript_id \"MUHGROUP2.$read\"\n";
			}
		}
		$read{$read} = "$begGTF$midGTF$endGTF";
	}
	close $in1;
	print "\n";
	print "Output = $LCY$input1.gtf$N\n";
	my $count = 0;
	my $groupcount =1;
	if (keys %{$dupe{grp}} == 0) {
		next;
	}
	my $totalcount = (keys %read);
	print $out1 "track type=gtf name='$label ($totalcount) $fileName1' color=$color visibility=0\n";
	foreach my $group (sort {$a <=> $b} keys %{$dupe{grp}}) {
		for (my $i = 0; $i < @{$dupe{grp}{$group}}; $i++) {
			my $read = $dupe{grp}{$group}[$i]; next if not defined $read{$read};
			$groupcount ++ if $i == 0;
			my $line = $read{$read};
			$line =~ s/MUHGROUP2/$groupcount/g;
#die "Undefined read=$read from group $group\n" if not 
			print $out1 "$line\n";
#			last if $count > 20;
			$count ++;
		}
#		last if $count > 20;
	}
}
}
__END__


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

# 0 is bad or not data (6)
# 1 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak
