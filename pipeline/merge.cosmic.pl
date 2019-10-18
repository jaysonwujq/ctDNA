use strict;
use FindBin qw($Bin $Script);
open COSMIC, "$Bin/../data_base/all.cosmic.tsv";
open RS,     "$Bin/../data_base/target.alf";
open IN,     "$ARGV[0]";
open OUT,    ">$ARGV[0].out";
my %cosmic;
my @as;
my $line;
my %phe_FATHMM_prediction;
my %rs;
my %rs_ratio;

while (<COSMIC>) {
  chomp;
  if ( $_ !~ "#" ) {
    @as = split;
    $cosmic{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] } = $as[2];
    $phe_FATHMM_prediction{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] }
      = "$as[5]\t$as[6]";
  }
}
my %real_snv;
while (<RS>) {
  chomp;
  @as                                               = split;
  $rs{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] }       = $as[2];
  $rs_ratio{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] } = $as[5];
}

while (<IN>) {
  chomp;
  if ( $_ =~ "#" ) {
    print OUT "$_\n";
  }
  else {
    @as = split /\t/, $_;
    if ( exists $cosmic{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] } ) {
      $as[5] = "$cosmic{$as[0]}{$as[1]}{$as[3]}{$as[4]}";
      $as[-1]
        = "$as[-1]\t$phe_FATHMM_prediction{$as[0]}{$as[1]}{$as[3]}{$as[4]}";
    }
    else {
      $as[-1] = "$as[-1]\tNA\tNA";
    }
    if ( exists $rs{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] } ) {
      $as[2]  = $rs{ $as[0] }{ $as[1] }{ $as[3] }{ $as[4] };
      $as[-1] = "$as[-1]\t$rs_ratio{$as[0]}{$as[1]}{$as[3]}{$as[4]}";
    }
    else {
      $as[-1] = "$as[-1]\tNA";
    }
    $line = join "\t", @as;
    print OUT "$line\n";
  }
}
