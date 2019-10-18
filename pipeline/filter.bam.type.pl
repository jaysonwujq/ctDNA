#!/bin/usr/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $indir, $sample, $target );
GetOptions(
  "in|i:s"     => \$indir,
  "sample|s:s" => \$sample,
  "target:s"   => \$target,
);
my $bin = $Bin;
if ( !$target){ $target = "$bin/../data_base/first.panel.bed";}
open my $DP,     "<", "$indir/bam/$sample.bam.filter";
open my $ANNO,   "<", "$indir/result/$sample.somtic.xls";
open my $OUT,    ">", "$indir/result/$sample.somtic.new.xls";
open my $LOSS,   ">", "$indir/result/$sample.somtic.mayloss.xls";
open my $NCBI,   ">", "$indir/result/$sample.somtic.ncbi.xls";
open my $TRANS,  "<", "$bin/../data_base/tran.id.new4";
open my $TARGET, "<", "$target";
open my $DPALL,  "<", "$indir/DP/$sample.target.bam.depth";
my %hash;
my %type;
my %map_pos;
my %loss;
my %testloss;
my %target;
my %dp;

while (<$TARGET>) {
  chomp;
  my ( $chr, $start, $end ) = split;
  foreach my $pos ( $start - 10 .. $end + 10 ) {
    $target{$chr}{$pos} = 1;
  }
}

while (<$DPALL>) {
  chomp;
  my ( $chr, $pos, $dp ) = split;
  $dp{$chr}{$pos} = $dp;
}

while (<$DP>) {
  chomp;
  my @dp = split /\t|\//, $_;
  if ( $dp[2] eq "-" )
  { ##insertion pos is not along to vcf file;because vcf is the insert pos before ,my pos is after
    $dp[1] = $dp[1] - 1;
  }
  $map_pos{ $dp[0] }{ $dp[1] }{ $dp[2] }{ $dp[3] } = join "\;",
    @dp[ 9 .. $#dp ];
  $hash{ $dp[0] }{ $dp[1] }{ $dp[2] }{ $dp[3] } = join "\t", @dp[ 4 .. 8 ];
  if ( ( $dp[5] >= 5 ) && ( $dp[6] >= 1 ) && ( $dp[7] >= 1 ) ) {
    my $key = join "\t", @dp[ 0 .. 3 ];
    $testloss{$key} = "$dp[4]\t$dp[5]";
    $type{ $dp[0] }{ $dp[1] }{ $dp[2] }{ $dp[3] } = $dp[5];

    #$map_pos{$dp[0]}{$dp[1]}{$dp[2]}{$dp[3]}=join "\;",@dp[9..$#dp];
  }
}
my %ncbi;
my %ensemble;
while (<$TRANS>) {
  chomp;
  my @line = split;
  foreach my $i ( 1 .. 2 ) {
    if ( !$line[$i] ) {
      $line[$i] = "NA";
    }
  }
  $ensemble{ $line[0] } = $line[1];
  $ncbi{ $line[0] }     = $line[2];
}
my @tmp;
my @anno;
print $OUT
  "Chr\tStart\tEnd\tRef\tAlt\trs\t1000_alt_ratio\tcosmic_id\tcosmic-FATHMM-prediction\tcosmic-FATHMM_prediction-score\tGene.refGene\tGeneDetail.ref\tGeneExonicFunc.refGene\t\t\tRef-dp\talt_dp\talt_ratio\talt_forward-dp\talt_reverse-dp\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tfilter-alt_dp\tfilter-alt_type\tfilter-alt_forward_type\tfilter-alt_reverse_type\tfilter-alt-both_forward_and_reverse\tfilter-alt_reads_map_pos\tis_filered\tTranscript-Ensemble-id\tgene\tTranscript-NCBI-id\ttexon\tcDNA-var\tAmino_acids_var\n";
while (<$ANNO>) {
  chomp;
  my @tmp = split /\t/, $_;
  my $all = join "\t", @tmp[ 0 .. 25 ];
  my $anno1 = "NA\tNA\tNA\tNA\tNA";
  if ( ( exists $hash{ $tmp[0] }{ $tmp[1] }{ $tmp[3] }{ $tmp[4] } )
    && ( exists $type{ $tmp[0] }{ $tmp[1] }{ $tmp[3] }{ $tmp[4] } ) ) {
    if ( $tmp[13] =~ "\:" ) {
      @anno = split /\:|\,/, $tmp[13];
      if ( !exists $ensemble{ $anno[0] } ) {
        $ensemble{ $anno[0] } = "NA";
      }
      for ( my $i = 0; $i <= $#anno; $i += 5 ) {
        if ( ( exists $ncbi{ $anno[$i] } )
          && ( $anno[ $i + 1 ] eq "$ncbi{$anno[$i]}" ) ) {
          $anno1 = join "\t", @anno[ $i .. $i + 4 ];
          print $NCBI "$anno1\n";
        }
      }
      my $anno = join "\t", @anno;
      print $OUT
        "$all\t$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\t$map_pos{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\tok\t$ensemble{$anno[0]}\t$anno1\t$anno\n";
      my $key = "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]";
      $loss{$key} = 1;
    }
    else {
      print $OUT
        "$all\t$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\t$map_pos{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\tok\n";
      my $key = "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]";
      $loss{$key} = 1;
    }
  }
  elsif  ( !exists $hash{ $tmp[0] }{ $tmp[1] }{ $tmp[3] }{ $tmp[4] } ){
    if  (( $tmp[3] =~ "-" ) || ( $tmp[4] =~ "-" ) )  {
      if ( $tmp[13] =~ "\:" ) {
        @anno = split /\:|\,/, $tmp[13];
        my $anno = join "\t", @anno;

        for ( my $i = 0; $i <= $#anno; $i += 5 ) {
          if ( ( exists $ncbi{ $anno[$i] } )
            && ( $anno[ $i + 1 ] eq "$ncbi{$anno[$i]}" ) ) {
            $anno1 = join "\t", @anno[ $i .. $i + 4 ];
          }
        }
        if ( !exists $ensemble{ $anno[0] } ) {
          $ensemble{ $anno[0] } = "NA";
        }
        print $OUT
          "$all\tNA\tNA\tNA\tNA\tNA\tNA\tnone-filter\t$ensemble{$anno[0]}\t$anno1\t$anno\n";
      }
      else {
        print $OUT "$all\tNA\tNA\tNA\tNA\tNA\tNA\tnone-filter\n";
      }
    }
    else{
      print $OUT "$all\tNA\tNA\tNA\tNA\tNA\tNA\tfiltered\n";
    }
  }
  elsif (exists $hash{ $tmp[0] }{ $tmp[1] }{ $tmp[3] }{ $tmp[4] } ) {

    if ( $tmp[13] =~ "\:" ) {
      @anno = split /\:|\,/, $tmp[13];
      if ( !exists $ensemble{ $anno[0] } ) {
        $ensemble{ $anno[0] } = "NA";
      }
      for ( my $i = 0; $i <= $#anno; $i += 5 ) {
        if ( ( exists $ncbi{ $anno[$i] } )
          && ( $anno[ $i + 1 ] eq "$ncbi{$anno[$i]}" ) ) {
          $anno1 = join "\t", @anno[ $i .. $i + 4 ];
        }
      }
      my $anno = join "\t", @anno;
      print $OUT
        "$all\t$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\t$map_pos{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\tfiltered\t$ensemble{$anno[0]}\t$anno1\t$anno\n";
      my $key = "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]";
      $loss{$key} = 1;
    }
    else {
      print $OUT
        "$all\t$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\t$map_pos{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}\tfiltered\n";
      my $key = "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]";
      $loss{$key} = 1;

    }
  }
}
foreach my $key ( keys %testloss ) {
  if ( !exists $loss{$key} ) {
    my @tmp  = split /\t/, $testloss{$key};
    my @keys = split /\t/, $key;
    if ( ( $tmp[0] > 15 )
      && ( $tmp[1] > 7 )
      && ( exists $target{ $keys[0] }{ $keys[1] } ) ) {
      print $LOSS
        "$key\t$hash{$keys[0]}{$keys[1]}{$keys[2]}{$keys[3]}\t$dp{$keys[0]}{$keys[1]}\n";
    }

    #    print "$key\t$hash{$key}\t$type{$key}\n";
  }
}
