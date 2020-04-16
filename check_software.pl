#!/usr/bin/perl
###

use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);

BEGIN {
	my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";
	my ($samtools) = `samtools 2>&1 | grep Version`; $samtools = "Unknown Samtools Version!" if not defined $samtools;
	my ($bedtools) = `bedtools --version`;  $bedtools = "Unknown bedtools Version!" if not defined $bedtools;
	my ($bowtie2) = `bowtie2 --version | grep version`;  $bowtie2 = "Unknown bowtie2 Version!" if not defined $bowtie2;
	my ($bismark) = `bismark --version| grep Version`;  $bismark = "Unknown bismark Version!" if not defined $bismark;
	my ($bismark_genome_preparation) = `bismark_genome_preparation --version | grep Version`;  $bismark_genome_preparation = "Unknown bismark_genome_preparation Version!" if not defined $bismark_genome_preparation;
	my ($R) = `R --version |grep version`; $R = "Unknown R Version!" if not defined $R;

	if ($samtools !~ /Version: (0\.1\.(19|[2-9]\d*)|0\.2|[1-9])/i) {
		print "Please install samtools at least version 0.1.19 before proceeding!\n\nsamtools=$samtools\n\n";
		$samtools = 0;
	}
	if ($bedtools !~ /bedtools v(2\.(2[5-9]|[3-9]\d*)|[3-9])/) {
		print "Please install bedtools at least version 2.25.0 before proceeding!\n\nbedtools=$bedtools\n\n";
		$bedtools = 0;
	}
	if ($bowtie2 !~ /version [2-9]\./) {
		print "Please install bowtie2 at least version 2.2.6 before proceeding!\n\nbowtie2=$bowtie2\n\n";
		$bowtie2 = 0;
	}
	if ($bismark !~ /v?(0\.2[0-9])/) {
		print "Please install bismark at least version 0.20.0 before proceeding!\n\nbismark=$bismark\n\n";
		$bismark = 0;
	}
	if ($bismark_genome_preparation !~ /v?(0\.[2-9])/) {
		print "\n\nPlease install bismark_genome_preparation at least version 0.20.0 before proceeding!\n\n\nbismark_genome_preparation=$bismark_genome_preparation\n\n";
		$bismark_genome_preparation = 0;
	}
	if ($R !~ /version (3\.(4\.[4-9]|[5-9])|[4-9])/) {
		print "Please install R at least version 3.4.4 before proceeding!\n\nR=$R\n\n";
		$R = 0;
	}
	my ($samtools_version) = $samtools eq 0 ? "NA" : $samtools =~ /Version.\s+(.+)$/; die "Can't determine samtools version from\n$samtools\n\n" if not defined $samtools_version;
	my ($bedtools_version) = $bedtools eq 0 ? "NA" : $bedtools =~ /bedtools\s+(.+)$/; die "Can't determine bedtools version from\n$bedtools\n\n" if not defined $bedtools_version;
	my ($bowtie2_version) = $bowtie2 eq 0 ? "NA" : $bowtie2 =~ /version\s+(.+)$/; die "Can't determine bowtie2 version from\n$bowtie2\n\n" if not defined $bowtie2_version;
	my ($bismark_version) = $bismark eq 0 ? "NA" : $bismark =~ /Version.\s+(.+)$/; die "Can't determine bismark version from\n$bismark\n\n" if not defined $bismark_version;
	my ($bismark_genome_preparation_version) = $bismark_genome_preparation eq 0 ? "NA" : $bismark_genome_preparation =~ /Version.\s+(.+)$/; die "Can't determine bismark_genome_preparation version from\n$bismark_genome_preparation\n\n" if not defined $bismark_genome_preparation_version;
	my ($R_version) = $R eq 0 ? "NA" : $R =~ /version\s+(\d+\.\d+\.\d+)\s*/; die "Can't determine R version from\n$R\n\n" if not defined $R_version;
	
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
	
	if ($samtools eq 0 or $bedtools eq 0 or $bowtie2 eq 0 or $bismark eq 0 or $bismark_genome_preparation eq 0 or $R eq 0) {
		print "\nPlease make sure to install these softwares first!\n";	
		print "bedtools (v2.25.0), bowtie2 (v2.2.6), bismark2 (v0.20.0), R (v3.4.4)\n\n";
	}

	print "\n--------------------\n\e[1;33m Software Check\e[0m\n--------------------\n\n";
	print "- samtools \e[1;36m$samtools_version\e[0m exists: \e[1;32m" . `which samtools` . "\e[0m" if $samtools ne 0;
	print "- bedtools \e[1;36m$bedtools_version\e[0m exists: \e[1;32m" . `which bedtools` . "\e[0m" if $bedtools ne 0;
	print "- bowtie2 \e[1;36m$bowtie2_version\e[0m exists: \e[1;32m" . `which bowtie2` . "\e[0m" if $bowtie2 ne 0;
	print "- bismark \e[1;36m$bismark_version\e[0m exists: \e[1;32m" . `which bismark` . "\e[0m" if $bismark ne 0;
	print "- bismark_genome_preparation \e[1;36m$bismark_genome_preparation_version\e[0m exists: \e[1;32m" . `which bismark_genome_preparation` . "\e[0m" if $bismark_genome_preparation ne 0;
	print "- R \e[1;36m$R_version\e[0m exists: \e[1;32m" . `which R` . "\e[0m" if $R ne 0;
	print "\n- Pushed \e[1;36m$libPath\e[0m into perl library paths at \@INC\n";
}

use myFootLib; use FAlite;

# MD5 script check
my $md5script = `which md5` =~ /md5/ ? "md5" : `which md5sum` =~ /md5sum/ ? "md5sum" : die "md5 or md5sum does not exist!\n";
print "- md5(sum) script: \e[1;36m$md5script\e[0m\n";
# Home Dir
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
die "Can't find $footLoopScriptsFolder\n" if not -d $footLoopScriptsFolder;
print "- FootLoop script folder: \e[1;36m$footLoopScriptsFolder\e[0m\n";
print "\n--------------------\n";

# Version
my @version = `cd $footLoopScriptsFolder && git log | head -n 15`;
my $version; my $commit = 0;
foreach my $line (@version[0..@version-1]) {
	last if $commit >= 2;
	chomp($line);
	if ($line =~ /^commit \w+$/) {
		$commit ++;
	}
	if ($line =~ /^\s+V\d+\.?\d*\w*\s*/i and $commit < 2) {
		my ($version2) = $line =~ /^\s+(V\d+\.?\d*\w*)\s*/i;
		$version .= "$version2\n";
	}
	elsif ($line !~ /^\s*$/ and $commit < 2) {
		$version .= "$line\n";
	}
}
$version = "UNKNOWN VERSION" if not defined $version;
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
	last if $version_small ne "vUNKNOWN";
}

print "\n--------------------\n\e[1;33m footLoop Version $version_small\e[0m\n--------------------\n\n";
print "$version\n---------------\n\n";

