use strict;
open my $LIST,   "<", "$ARGV[0]";
open my $TARGET, "<", "$ARGV[1]";
open my $OUT,    ">", "$ARGV[2]";
my $in = $ARGV[3];
my %hash;
my @as;
my $i = 0;
while (<$TARGET>) {
  chomp;
  @as = split;
  foreach ( $as[1] .. $as[2] ) {
    $hash{ $as[0] }{$_} = 1;
  }
}
my $num;
my $sum;
my ($dp2000, $dp200, $dp100,$dp500);
my $sum_rmdup;
my $num_rmdup;
print $OUT "Sample\tdp_mean_raw\tdp>100.ratio\tdp>200.ratio\tdp>500.ratio\tdp>2000.ratio\tdp_mean_rmdup\n";
my $count=0;
while (<$LIST>) {
  chomp;
  $count++;
  $num       = 0;
  $sum       = 0;
  $dp2000    = 0;
  $dp200     = 0;
  $dp100     = 0;
  $dp500     = 0;
  $sum_rmdup = 0;
  $num_rmdup = 0;
  my $sample = $_;
  open my $DP_SAM, "<", "$in/$sample.target.bam.depth";

  while (<$DP_SAM>) {
    @as = split;
    if ( exists $hash{ $as[0] }{ $as[1] } ) {
      $num++;
      $sum += $as[2];
      if ( $as[2] > 2000 ) {
        $dp2000++;
      }
      if ( $as[2] > 200 ) {
        $dp200++;
      }
      if ( $as[2] > 100 ) {
        $dp100++;
      }
      if ( $as[2] > 500 ) {
        $dp500++;
      }
    }
  }
  print "$count\t$sum\t$num\n";
  my $dp_raw        = $sum / $num;
  my $dp_2000_ratio = $dp2000 / $num;
  my $dp_200_ratio = $dp200 / $num;
  my $dp_100_ratio = $dp100 / $num;
  my $dp_500_ratio = $dp500 / $num;
  print $OUT "$sample\t$dp_raw\t$dp_100_ratio\t$dp_200_ratio\t$dp_500_ratio\t$dp_2000_ratio\n";
}
