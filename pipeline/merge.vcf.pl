#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out );
GetOptions(
  "in:s"  => \$in,
  "out:s" => \$out,
);

open my $LIST, "<", "$in";
open my $VCF,  ">", "$out";
my $i = 0;
while (<$LIST>) {
  chomp;
  $i++;
  open my $IN, "<", "$_";
  while (<$IN>) {
    chomp;
    if ( $_ =~ /^#/ ) {
      if ( $i == 1 ) {
        print $VCF "$_\n";
      }
      else {
        next;
      }
    }
    else {
      print $VCF "$_\n";
    }
  }
}
