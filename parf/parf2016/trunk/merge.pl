#!/usr/bin/perl
use File::Basename;
use POSIX "floor";

$first = shift;
$second = shift;

unless ($second) {
print <<EOF;
Usage: splitset.pl <forest1> <forest2>

This utility merges forest2 into forest1. 
EOF
  exit 0;
}

sub log10 {
  my $n = shift;
  return log($n) / log(10);
}

sub exp10 {
  my $n = shift;
  return exp($n * log(10));
}

sub pad {
  my $n = shift;
  my $l = shift;
  return sprintf("%0".$l."d", $n);
}

$first = $1 if ($first =~ /^(.*)\.forest/);
$second = $1 if ($second =~ /^(.*)\.forest/);

$firstdir = dirname($first);
$seconddir = dirname($second);

open(F, "<$second.forest") || die("Error: Can't open $second.forest");
$numsecond = <F>;
close(F);

rename("$first.forest", "$first.forest.tmp");
open(F, "<$first.forest.tmp") || die("Error: Can't open $first.forest.tmp");
open(R, ">$first.forest") || die("Error: Can't open $first.forest");
while (<F>) {
  if (!$numfirst) {
    $numfirst = $_;
    print R ($numfirst + $numsecond), "\n";
  } else {
    print R;
  }
}
close(R);
close(F);
unlink("$second.forest");
unlink("$first.forest.tmp");

$sizefirst = floor(log10($numfirst)) + 1;
$sizesecond = floor(log10($numsecond)) + 1;
$sizeresult = floor(log10($numfirst + $numsecond)) + 1;

$limit = exp10($sizeresult - 1) - 1;
for ($i = 1; $i <= $limit; $i++) {
  rename("$first.".pad($i, $sizefirst).".tree",
    "$first.".pad($i, $sizeresult).".tree");
}

for ($i = 1; $i <= $numsecond; $i++) {
  rename("$second.".pad($i, $sizesecond).".tree",
    "$first.".pad($i+$numfirst, $sizeresult).".tree");
}
