#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $chr );
use List::Util qw(max);
use List::Util qw(sum);
GetOptions(
  "in|i:s"  => \$in,
  "out|o:s" => \$out,
);

open my $IN,  "<", "$in";
open my $OUT, ">", "$out";

my $key    = "";
my $before = "";
my @stats;
my $line;
my $each_num = 18;
while ( $line = <$IN> ) {
  chomp $line;
  my @tmp = split /\t/, $line;
  $key = join "\t", @tmp[ 0 .. 4 ];
  if ( $key eq $before ) {
    push @stats, $line;
  }
  else {
    my $tmp = join "\t", @stats;
    my $result = &rmdup($tmp);
    print $OUT "$result\n";
    @stats = split /\t/, $line;
    $before = $key;
  }
}

sub rmdup {
  my ($tmp) = $_[0];
  my @tmp     = split /\s+/, $tmp;
  my $num     = ( $#tmp + 1 ) / $each_num;
  my $norm1   = 0;
  my $alt1    = 0;
  my $norm2   = 0;
  my $alt2    = 0;
  my $quality = 0;
  my %alt;
  my $reads_name;

  foreach my $i ( 0 .. $num - 1 ) {
    if (
      ( $tmp[ $i * $each_num + 5 ] + $tmp[ $i * $each_num + 6 ] ) > $quality )
    {
      $reads_name = $tmp[ $i * $each_num + 11 ];
      $quality    = $tmp[ $i * $each_num + 5 ] + $tmp[ $i * $each_num + 6 ];
    }
  }
  return $reads_name;
}
