#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_n);
getopts("vn:");

my ($footPeak_folder) = ($opt_n);
die "\nusage: $YW$0$N -n $CY<footPeak folder>$N\n\n" unless defined $opt_n;

my @folders = split("/", $footPeak_folder);
my $outName = $folders[@folders-1]; $outName =~ s/\///g;

my $cmd = "cat $footPeak_folder/99_FOOTSTATS/.PEAKSTATSTEMP/.0_RESUL*TXT|perl -pi -e 's/[ ]+/_/g' |grep -vP '^#' |perl -pi -e 's/^.+gene(.+)_(Pos|Neg)_.+\\t(.+)\\t([A-Z][A-Z])\\t(\\d+)\\t(.+)\$/\$1\\t\$2\\t\$4\\t\$5\\t\$6/' > $outName.TEMP";
system("$cmd") == 0 or die "\n\nFailed $LCY$cmd$N: $!\n\n";

print "\n$LCY$outName.TEMP$N\n\n";

__END__
cat PBEH2_240401_p1000_q0_L1_t40_l50_w15/99_FOOTSTATS/.PEAKSTATSTEMP/.0_RESUL*TXT|perl -pi -e 's/[ ]+/_/g' |grep -vP "^#" |perl -pi -e 's/^.+gene(.+)_(Pos|Neg)_.+\t(.+)\t([A-Z][A-Z])\t(\d+)\t(.+)$/$1\t$2\t$4\t$5\t$6/' > 
CELL2_240325_240401_t40w15l50.TEMP

