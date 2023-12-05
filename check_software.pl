#!/usr/bin/perl
###

use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);

BEGIN {
	my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
	push(@INC, $libPath);
}

use myFootLib; use FAlite;
my ($nocolor) = @ARGV;
my $warningslog = "";
my $softwareslog = "";
my $softwaresloglong = "";

my @softwares = qw(bismark bedtools bowtie2 samtools Rscript);
foreach my $software (@softwares[0..@softwares-1]) {
	my ($location) = `which $software`;
	my $versionprint = `$software --version 2>&1`;
	if (not defined $location or (defined $location and $location eq "")) {
		$warningslog .= "$software doesn't exist!\n";
		$location = "N/A";
		$versionprint="N/A";
		$softwareslog .= "$software=$location\n";

	}
	else {
		chomp($location);
		chomp($versionprint);	
		$softwareslog     .= "$software=$location\n";
	}
	$softwaresloglong .= "#$software --version\n$LGN$versionprint\n\n";
}
# MD5 script check
my $md5script = `which md5`    =~ /md5/    ? `which md5`    : 
				 	 `which md5sum` =~ /md5sum/ ? `which md5sum` :
					 die "md5 or md5sum does not exist!\n";
chomp($md5script);

# Home Dir
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
die "Can't find footLoop folder $LCY$footLoopScriptsFolder\n" if not -d $footLoopScriptsFolder;

# FootLoop Version
my @version = `cd $footLoopScriptsFolder && git log`;
my ($version, $commit, $commitdate, $versiondate, $commitnum);
$commitnum = 0;
my $linecount = 0;
foreach my $line (@version[0..@version-1]) {
	chomp($line);
	$linecount ++;
	if ($line =~ /commit .+$/) {
		$commitnum ++;
		next if defined $commit;
		($commit) = $line =~ /commit (.+)$/;
	}
	elsif ($line =~ /^\s+[Vv]\d+\.?\d*\w*\s*/i) {
		($version) = $line =~ /^\s+([Vv]\d+\.?\d*\w*)\s*/i;
	}
	elsif ($line =~ /^Date: /) {
		($commitdate) = $line =~ /^Date:\s+(\w+.+)$/ if not defined $commitdate;
		($versiondate) = $line =~ /^Date:\s+(\w+.+)$/;
	}
	last if defined $version;
}
$versiondate =~ s/[ :]/_/g; $versiondate =~ s/\-(\d+)/min$1/;
$commitdate =~ s/[ :]/_/g; $commitdate =~ s/\-(\d+)/min$1/;

my $scripttestlog = "$version\_$commit
footLoop_version=$version\_$commit
version=$version
commit=$commit
commitdate=$commitdate
versiondate=$versiondate
footLoop_script_folder=$footLoopScriptsFolder
md5sum_script=$md5script

$YW# Locations:
$softwareslog

$YW# Versions:
$softwaresloglong

$version\_$commit\n";
$warningslog = "${LRD}\nWARNING!!:\n$warningslog\n$N" if $warningslog ne "";
print "$scripttestlog";
print STDERR "$warningslog\n";
__END__
#my $bismarkloc  = `which bismark` ; chomp($bismarkloc) ; $bismarkloc  = "N/A" if $bismarkloc eq "";
#my $bedtoolsloc = `which bedtools`; chomp($bedtoolsloc); $bedtoolsloc = "N/A" if $bedtoolsloc eq "";
#my $bowtie2loc  = `which bowtie2` ; chomp($bowtie2loc) ; $bowtie2loc  = "N/A" if $bowtie2loc eq "";
#my $samtoolsloc = `which samtools`; chomp($samtoolsloc); $samtoolsloc = "N/A" if $samtoolsloc eq "";
#my $Rscriptloc  = `which Rscript` ; chomp($Rscriptloc) ; $Rscriptloc  = "N/A" if $Rscriptloc eq "";
#my $bismarktest  = `bismark  --version 2>&1 | head`; chomp($bismarktest);
#my $bedtoolstest = `bedtools --version 2>&1 | head`; chomp($bedtoolstest);
#my $bowtie2test  = `bowtie2  --version 2>&1 | head`; chomp($bowtie2test);
#my $samtoolstest = `samtools --version 2>&1 | head`; chomp($samtoolstest);
#my $Rscripttest  = `Rscript  --version 2>&1 | head`; chomp($Rscripttest);


#$version = "UNKNOWN VERSION" if not defined $version;
#my ($version_small) = "vUNKNOWN";
#foreach my $versionz (@version[0..@version-1]) {
#   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
#	last if $version_small ne "vUNKNOWN";
#}

bedtools=$bedtoolsloc
samtools=$samtoolsloc
bowtie2=$bowtie2loc
bismark=$bismarkloc
Rscript=$Rscriptloc
#bedtools --version
$LGN$bedtoolstest

#samtools --version
$LGN$samtoolstest

#bowtie2 --version
$LGN$bowtie2test

#bismark --version
$LGN$bismarktest

#Rscript --version
$LGN$Rscripttest
