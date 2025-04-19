#!/usr/bin/perl
###

use warnings; use strict; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . 'lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";
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
