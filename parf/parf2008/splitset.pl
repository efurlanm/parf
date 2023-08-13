#!/usr/bin/perl

$set = shift;
$numtest = shift;

unless ($set) {
print <<EOF;
Usage: splitset.pl <set> [<ntest>]

If called with two parameters, splitset.pl splits the
instance set into a trainset of (size - ntest) instances
and a testset of (ntest) instances.

If called with only one parameter, splitset.pl will count
the number of instances in the set.
EOF
  exit 0;
}
$set =~ s/\.arff$//;

open(R, "<$set.arff");
$data = 0;
while (<R>) {
  $numall++ if ($data);
  next if (/^\s*(%|$)/);
  $data = 1 if (/\@data/);
}
close(R);

if ($numtest) {
  open(R, "<$set.arff");
  open(TR, ">$set-train.arff");
  open(TS, ">$set-test.arff");

  $data = 0;
  for ($i = 0; $i < $numtest; $i++) {
    $n = int(rand($numall)) while ($test[$n]);
    $test[$n] = 1;
  }
  $i = 0;
  $data = 0;
  while (<R>) {
    if ($data) {
      next if (/^\s*(%|$)/);
      if ($test[$i++]) {
        print TS;
      } else {
        print TR;
      }
    } else {
      print TR;
      print TS;
    }
    $data = 1 if (/\@data/);
  }
  $numtrain = $numall - $numtest;
  print "Split $numall into $numtrain training and $numtest test\n";
} else {
  print "Total instances: $numall\n";
}
