use strict;

open my $LIST, "<", "$ARGV[0]";
my $in = $ARGV[1];
open my $OUT, ">>", "$ARGV[2]";

print $OUT "sample\tgene\tpos_num\tlog2_ratio\tcopy_number\n";
while (<$LIST>) {
  chomp;
  my $name = $_;
  my @tmp  = split /\./,$name;
  my $sample = $tmp[0];
  my $gene   = $tmp[1];
  open my $IN, "<", "$in/$name";
  my $all     = 0;
  my $cnv     = 0;
  my $cnv_sum = 0;
  while (<$IN>) {
    chomp;
    $all++;
    my @as = split /\s+/, $_;
    if ( $as[1] > 0.7 ) {
      $cnv++;
      $cnv_sum += $as[1];
    }
  }
  my $ratio    = 0;
  my $cnv_num  = 0;
  my $copy_num = 0;
  if ( $cnv > 0 ) {
    $ratio    = $cnv / $all;
    $cnv_num  = $cnv_sum / $all;
    $copy_num = 2**$cnv_num+2;
  }
  if ( ( $ratio > 0.7 ) && ( $all > 500 ) ) {
    print $OUT "$sample\t$gene\t$all\t$cnv_num\t$copy_num\n";
  }
}
