#!/usr/bin/perl

$r = shift;
$w = shift;
$wp = $w;
$w = "$r.tmp" unless ($wp);
if ($r) {
  open (R, "<$r") or die "Error: Can't open $r for input\n";
  open (W, ">$w") or die "Error: Can't open $w for output\n";
} else {
  open (R, "<&STDIN") or die "Error: Can't open STDIN\n";
  open (W, ">&STDOUT") or die "Error: Can't open STDOUT\n";
}
while (<R>) {
  chomp;
  s/\cM$//g; # for DOS-encoded files (CRLF)
  s/\\'/"/g; # for weka quote escapes
  while (length($_) >= 1022 && ($p = rindex($_, ',', 1022)) != -1) {
    print W substr($_, 0, ++$p)."&\n";
    $_ = substr($_, $p);
  }
  print W "$_\n";
}
close W;
close R;
rename ($w, $r) if ($r && !$wp);

