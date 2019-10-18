#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use warnings;

my ( $in, $out, $chr, $in1 );
use List::Util qw(max);
GetOptions(
  "in|i:s"  => \$in,
  "bam:s"   => \$in1,
  "out|o:s" => \$out,
);
use Bio::DB::Sam;
my $sam = Bio::DB::Sam->new(
  -bam   => "$in",
  -fasta => "/data1/software/b37/human_g1k_v37.fasta",
);
open my $OUT,  ">", "$out";
open my $LOG,  ">", "$out.filter";
open my $LOG1, ">", "$out.log";
my (
  $pair,       $dup1,         $dup2,    $nm1,     $nm2,       $md1,
  $first_mate, $second_mate,  $start1,  $start2,  $seq_name,  $md2
);
my @pairs = $sam->get_features_by_location(
  -type => 'read_pair',
);
my $num = 0;
my %rmxa;
foreach my $pair (@pairs) {
  my %reads;
  ( $first_mate, $second_mate )
    = $pair->get_SeqFeatures;    #distinguish first and second reads

  if ( $first_mate && $second_mate ) {
    $start1 = $first_mate->start;     #first reads start mapping pos
    $start2 = $second_mate->start;    #second reads start mapping pos
    if ( ($start1) && ($start2) && ( $start1 > 0 ) && ( $start2 > 0 ) ) {
      $dup1 = $first_mate->get_tag_values('XA');   #first reads dup
      $dup2 = $second_mate->get_tag_values('XA');  #second reads dup
      $nm1  = $first_mate->get_tag_values('NM');   #First reads mismatch base numbers
      $nm2  = $first_mate->get_tag_values('NM');   #First reads mismatch base numbers
      $md1  = $first_mate->cigar_str;   #first reads mismatch base
      $md2  = $second_mate->cigar_str;  #second reads mismatch base
      $seq_name = $first_mate->name;
      my $xa = &filter_multialign($dup1, $dup2, $nm1, $nm2, $md1, $md2);
      if ($xa == 0 ){
        $rmxa{$seq_name} = 1;
      }
    }
  }
}

open my $BAM, "<", "$in1";
while (<$BAM>) {
  chomp;
  if ( $_ =~ /^@/ ) {
    print $OUT "$_\n";
  }
  else {
    my @tmp = split;
    my $key = $tmp[0];
    if (  !exists $rmxa{$key} )  {
      print $OUT "$_\n";
    }
  }
}

sub filter_multialign {
  my ( $xa1, $xa2, $nm1, $nm2, $md1, $md2 ) = @_;
  if ( $xa1 && $xa2 ) {
    my %class;
    my $xa1_nm1 = (split /\,|\;/, $xa1)[3];
    my $xa2_nm2 = (split /\,|\;/, $xa2)[3];
    my @nms     = ($nm1, $nm2, $xa1_nm1, $xa2_nm2);
    my @mds     = ($md1, $md2, $xa1, $xa2);
    my @names   = ("nm1", "nm2", "xa1_nm1", "xa2_nm2");
    foreach my $i (0..$#nms){
      if ($nms[$i] == 0 ){
        $class{$names[$i]} = 1;
      }elsif ($nms[$i] <= 2 ){
        $class{$names[$i]} = 2;
      }else{
        $class{$names[$i]} = 3;
      }
      if ($mds[$i] =~ /S/){
        $class{$names[$i]} ++;
      }
    }
    if( ($class{nm1} < $class{xa1_nm1}) || ($class{nm2} < $class{xa2_nm2}) ){
      print $LOG "non_filter_multialign\t$xa1, $xa2, $nm1, $nm2, $md1, $md2\n";
      return "1";
    }
    else{
      print $LOG "filter_multialign\t$xa1, $xa2, $nm1, $nm2, $md1, $md2\n";
      return "0";
    }
  }
  else {
    return "1";
  }
}
