#!/usr/bin/perl

use strict;
use warnings;

my ($string) = @ARGV;
die "Usage: $0 <stirng>\n" if @ARGV == 0;

my $length = length($string);
my $start = 0;
my $end = 2391;

while ($start < $length) {
  my $substring = substr($string, $start, $end);
  print "$substring\n";
  $start += 2392;
  $end += 2392;
  if ($end >= $length) {
    $end = $length - 1;
  }
}
