#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $list );
GetOptions(
  "in:s"  => \$in,
  "out:s" => \$out,
);

open my $IN,  "<", "$in/stats/list";
open my $OUT, ">", "$out/type";
open my $ALL, ">", "$out/all_types";

while (<$IN>) {
  chomp;
  my @tmp   = split /\_/, $_;
  my @types = split //,   $tmp[1];
  my $sample = $tmp[0];
  print $OUT "$types[2]\t$sample\t$types[0]\n";
  my $all = join "\t", @types;
  print $ALL "$sample\t$all\n";
}
