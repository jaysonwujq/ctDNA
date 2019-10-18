use strict;
open CLEAN,  "$ARGV[0]";
open BAMTAR, "$ARGV[1]";
open DP,     "$ARGV[2]";
open OUT,    ">$ARGV[3]";
my $i = 0;
while (<CLEAN>) {
  $i++;
  chomp;
  my $target = <BAMTAR>;
  my $dp     = <DP>;
  chomp $target;
  chomp $dp;
  my @clean  = split /\t/, $_;
  my @target = split /\t/, $target;
  my @dp     = split /\t/, $dp;

  if ( $i > 1 ) {
    $target[0] = $target[1] / $clean[3];
  }
  else {
    $target[0] = "mapping_rate";
  }
  $target = join "\t", @target;
  $dp     = join "\t", @dp[ 1 .. $#dp ];
  print OUT "$_\t$target\t$dp\n";
}
