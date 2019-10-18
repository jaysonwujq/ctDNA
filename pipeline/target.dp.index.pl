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
my $num2;
my $sum_rmdup;
my $num_rmdup;
print $OUT "Sample\tdp_mean_raw\tdp>2000.ratio\tdp_mean_rmdup\n";
while (<$LIST>) {
  chomp;
  $num       = 0;
  $sum       = 0;
  $num2      = 0;
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
        $num2++;
      }
    }
  }
  print "$sum\t$num\n";
  my @dps = (1,2,3,4);
  my @sum_rmdup;
  my @num_rmdup;
  foreach my $dp (@dps){
    open my $DP_RMDUP, "<", "$in/$sample.target.realigned.bam.$dp.depth";
    while (<$DP_RMDUP>) {
      @as = split;
      if ( exists $hash{ $as[0] }{ $as[1] } ) {
        $num_rmdup[$dp - 1] ++;
        $sum_rmdup[$dp - 1] += $as[2];
      }
    }
  }
#  open my $DP_RMDUP, "<", "$in/$sample.target.realigned.bam.depth";
#  print "$in/$sample.target.realigned.bam.depth\n";
#  while (<$DP_RMDUP>) {
#    @as = split;
#    if ( exists $hash{ $as[0] }{ $as[1] } ) {
#      $num_rmdup++;
#      $sum_rmdup += $as[2];
#    }
#  }
  print "$sum_rmdup\t$num_rmdup\n";
  my $dp_raw        = $sum / $num;
  my $dp_2000_ratio = $num2 / $num;
  #my $dp_rmdup      = $sum_rmdup / $num_rmdup;
  #print $OUT "$sample\t$dp_raw\t$dp_2000_ratio\t$dp_rmdup\n";
  print $OUT "$sample\t$dp_raw\t$dp_2000_ratio";
  foreach (0..3){
    my $dp_rmdup      = $sum_rmdup[$_] / $num_rmdup[$_];
    print $OUT "\t$dp_rmdup";
  }
  print $OUT "\n";
}
