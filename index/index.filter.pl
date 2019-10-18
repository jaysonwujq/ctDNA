#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $reads_num );
use List::Util qw(max);
use List::Util qw(sum);
GetOptions(
  "in|i:s"  => \$in,
  "out|o:s" => \$out,
  "num|n:s" => \$reads_num,
);
if (!$reads_num){$reads_num =1;}
open my $IN,     "<", "$in";
open my $OUT,    ">", "$out";
open my $LOG,    ">", "$out.log";
open my $FILTER, ">", "$out.filter";
my $key    = "";
my $before = "";
my @stats;
my $line;
my $each_num                 = 18;
my $all_reads                = 0;
my $contain_read             = 0;
my $all_type                 = 0;
my $contain_type             = 0;
my $only_one_read            = 0;
my $only_two_reads_diff_type = 0;
my %pos;
my %index;
my $may_wrong_index = 0;
my %may_wrong_index;

## count pos_types and index_types ;store the may wrong index (which only one read inport index but the pos have more than five reads)
while ( $line = <$IN> ) {
  chomp $line;
  my @tmp = split /\t/, $line;
  my $key_pos   = join "\t", @tmp[ 2 .. 4 ];
  my $key_index = join "\t", @tmp[ 0 .. 4 ];
  $pos{$key_pos}++;
  $index{$key_index}++;
}
close $IN;
my $only_one_pos   = 0;
my $all_pos_type   = 0;
my $all_index_type = 0;
foreach my $key_index ( keys %index ) {
  my @keys_all = split /\t/, $key_index;
  my $key_pos = join "\t", @keys_all[ 2 .. 4 ];
  if ( ( $index{$key_index} == 1 ) && ( $pos{$key_pos} > 5 ) ) {
    $may_wrong_index++;
    $may_wrong_index{$key_index} = 1;
  }
  $all_index_type++;
}
foreach my $key_pos ( keys %pos ) {
  if ( $pos{$key_pos} == 1 ) {
    $only_one_pos++;
  }
  $all_pos_type++;
}

open $IN, "<", "$in";
while ( $line = <$IN> ) {
  $all_reads++;
  chomp $line;
  my @tmp = split /\t/, $line;
  $key = join "\t", @tmp[ 0 .. 4 ];
  if ( exists $may_wrong_index{$key} ) {
    next;
  }
  if ( $key eq $before ) {
    push @stats, $line;
  }
  else {
    $all_type++;
    my $all_same_reads_merge = join "\t", @stats;
    my $result = &rmdup($all_same_reads_merge);
    if ( $result ne "0" ) {
      print $OUT "$result\n";
      $contain_type++;
    }
    @stats = split /\t/, $line;
    $before = $key;
  }
}

## count stats
my $contain_read_ratio             = $contain_read / $all_reads;
my $contain_type_ratio             = $contain_type / $all_type;
my $only_one_read_type_ratio       = $only_one_read / $contain_type;
my $only_two_reads_diff_type_ratio = $only_two_reads_diff_type / $all_type;
my $only_one_pos_ratio             = $only_one_pos / $only_one_read;
my $may_wrong_index_ratio          = $may_wrong_index / $only_one_read;
my $index_ratio                    = $all_index_type / $all_pos_type;
print $LOG
  "all_reads\tfiltered_read\tall_type\tcontain_type\tonly_one_read\tfiltered_read_ratio\tcontain_type_ratio\tonly_one_read_type_ratio\tonly_two_reads_diff_type_ratio\tonly_one_pos_ratio\tmay_wrong_index_ratio\tall_index_type\tall_pos_type\tindex_ratio\n";
print $LOG
  "$all_reads\t$contain_read\t$all_type\t$contain_type\t$only_one_read\t$contain_read_ratio\t$contain_type_ratio\t$only_one_read_type_ratio\t$only_two_reads_diff_type_ratio\t$only_one_pos_ratio\t$may_wrong_index_ratio\t$all_index_type\t$all_pos_type\t$index_ratio\n";

sub rmdup {
  my ($all_same_reads_merge) = $_[0];
  my @all_same_reads_merge_split = split /\t/, $all_same_reads_merge;
  my $num_reads = ( $#all_same_reads_merge_split + 1 ) / $each_num;
  my $quality   = 0;
  my %alt;
  my %quality;
  my %num;
  my %alt1;
  my %alt2;
  my $key_alt;
  my $reads_name;
  my $alt_reads1 = 0;
  my $alt_reads2 = 0;
  my @alt_read1;
  my @alt_read2;

  foreach my $i ( 0 .. $num_reads - 1 ) {
    $alt_reads1 = 0;
    $alt_reads2 = 0;
    ## for the first reads var type
    if ( $all_same_reads_merge_split[ $i * $each_num + 9 ] =~ /A|T|C|G/ ) {
      @alt_read1 = split /\,/,
        $all_same_reads_merge_split[ $i * $each_num + 9 ];
      my $alt_merge = 0;
      my $pos_s     = 0;
      if ( $alt_read1[1] eq "S" ) {
        $pos_s = $alt_read1[0];
      }
      for ( my $j = 1; $j <= $#alt_read1; $j += 2 ) {
        if ( $alt_read1[$j] =~ /\// ) {
          my $pos = $alt_read1[ $j - 1 ]
            - $pos_s + $all_same_reads_merge_split[ $i * $each_num + 7 ];
          my $ref_alt_base = $alt_read1[$j];
          $alt_merge = "$alt_merge\t$pos\t$ref_alt_base";
        }
      }
      $alt_reads1 = $alt_merge;
      $alt1{$alt_reads1}++;
    }
    ## for the second  reads var type
    if ( $all_same_reads_merge_split[ $i * $each_num + 10 ] =~ /A|T|C|G/ ) {
      @alt_read2 = split /\,/,
        $all_same_reads_merge_split[ $i * $each_num + 10 ];
      my $alt_merge = 0;
      my $pos_s     = 0;
      if ( $alt_read2[1] eq "S" ) {
        $pos_s = $alt_read2[0];
      }
      for ( my $j = 1; $j <= $#alt_read2; $j += 2 ) {
        if ( $alt_read2[$j] =~ /\// ) {
          my $pos = $alt_read2[ $j - 1 ]
            - $pos_s + $all_same_reads_merge_split[ $i * $each_num + 8 ];
          my $ref_alt_base = $alt_read2[$j];
          $alt_merge = "$alt_merge\t$pos\t$ref_alt_base";
        }
      }
      $alt_reads2 = $alt_merge;
      $alt2{$alt_reads2}++;
    }
    $key_alt = "$alt_reads1\n$alt_reads2";
    ## each var type save the highest base quality
    if ( exists $alt{$key_alt} ) {
      if (
        (
            $all_same_reads_merge_split[ $i * $each_num + 5 ]
          + $all_same_reads_merge_split[ $i * $each_num + 6 ]
        ) > $quality{$key_alt}
        ) {
        $alt{$key_alt} = $all_same_reads_merge_split[ $i * $each_num + 11 ];
        $quality{$key_alt} = $all_same_reads_merge_split[ $i * $each_num + 5 ]
          + $all_same_reads_merge_split[ $i * $each_num + 6 ];
      }
      $num{$key_alt}++;
    }
    else {
      $alt{$key_alt}     = $all_same_reads_merge_split[ $i * $each_num + 11 ];
      $quality{$key_alt} = $all_same_reads_merge_split[ $i * $each_num + 5 ]
        + $all_same_reads_merge_split[ $i * $each_num + 6 ];
      $num{$key_alt} = 1;
    }
  }
  my $max_ratio        = 0;
  my $second_max_ratio = 0;
  foreach my $keys ( keys %alt ) {    ## save the type ratio more than 50%
    if ( ( $num{$keys} >= $reads_num ) && ( $num{$keys} / $num_reads > 0.5 )) {
      my $ratio = $num{$keys} / $num_reads;
      if ( $num_reads == 1 ) {
        $only_one_read++;
      }
      return "$alt{$keys}\t$num_reads\t$ratio";
    }
    elsif ( $num{$keys} / $num_reads > $max_ratio ) {
      $second_max_ratio = $max_ratio;
      $max_ratio        = $num{$keys} / $num_reads;
    }
    elsif ( $num{$keys} / $num_reads > $second_max_ratio ) {
      $second_max_ratio = $num{$keys} / $num_reads;
    }
  }
  my $max_alt_reads1;
  my $max_alt_reads2;
  ## save the type ratio for each reads more than 50%
  foreach my $key_reads1 ( keys %alt1 ) {
    if (  ( $alt1{$key_reads1} >= $reads_num) && ( $alt1{$key_reads1} / $num_reads > 0.5 ) ){
      $max_alt_reads1 = $key_reads1;
    }
  }
  foreach my $key_reads2 ( keys %alt2 ) {
    if (( $alt2{$key_reads2} / $num_reads > 0.5 ) && ($alt2{$key_reads2} >= $reads_num) ){
      $max_alt_reads2 = $key_reads2;
    }
  }
  if ( ($max_alt_reads1) && ($max_alt_reads2) ) {
    my $key   = "$max_alt_reads1\n$max_alt_reads2";
    my $ratio = $num{$key} / $num_reads;
    return "$alt{$key}\t$num_reads\t$ratio";
  }
  ## save the type ratio more than 30% and the second max is 20% samller
  foreach my $keys ( keys %alt ) {
    if ( ( $num{$keys} / $num_reads >= 0.3 )
      && ( $num_reads >= 3 )
      && ( $num{$keys} / $num_reads == $max_ratio )
      && ( $max_ratio - $second_max_ratio > 0.2 ) 
      && ( $num{$keys} >= $reads_num)
    ) {
      my $ratio = $num{$keys} / $num_reads;
      return "$alt{$keys}\t$num_reads\t$ratio";
    }
  }
  if ( $num_reads == 2 ) {
    $only_two_reads_diff_type++;
  }
  $contain_read += $num_reads;
  foreach my $i ( 0 .. $num_reads - 1 ) {    ## print the filtered type reads
    my $result = join "\t",
      @all_same_reads_merge_split[ ( $i * $each_num )
      .. ( $i * $each_num + $each_num - 1 ) ];
    print $FILTER "$result\n";
  }
  print $FILTER "\n";
  return "0";
}
