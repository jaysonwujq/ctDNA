#!/usr/bin/perl
use utf8;
open LIST, "$ARGV[0]";

open IN,  "gzip -dc /data1/workdir/liubei/program/cds/rs.POS.ALL.gz|";
open OUT, ">$ARGV[0].RS";
my %hash;
while (<LIST>) {
  chomp;
  my @as = split /\s+/, $_;
  $hash{ $as[2] } = "NA\tNA\tNA\tNA";
}
close LIST;
open LIST, "$ARGV[0]";
while (<IN>) {
  chomp;
  my @as = split;
  if ( exists $hash{ $as[0] } ) {
    $hash{ $as[0] } = "$as[1]\t$as[2]\t$as[3]\t$as[4]";
  }
}
while (<LIST>) {
  chomp;
  my @as = split /\s+/, $_;
  if ( exists $hash{ $as[2] } ) {
    print OUT "$_\t$hash{$as[2]}\n";
  }
}
