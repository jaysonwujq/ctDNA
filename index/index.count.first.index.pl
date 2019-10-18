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
  "chr:s"   => \$chr
);
use Bio::DB::Sam;
my $sam = Bio::DB::Sam->new(
  -bam   => "$in",
  -fasta => "/data1/software/b37/human_g1k_v37.fasta",
);
open my $OUT,  ">", "$out.chr$chr";
open my $LOG,  ">", "$out.chr$chr.filter";
open my $LOG1, ">", "$out.chr$chr.log";
my (
  $pair,    $length,      $seqid1,      $f_start, $start1,
  $start2,  $reads1,      $reads2,      $strand1, $strand2,
  $paired1, $paired2,     $md1,         $md2,     $dup1,
  $dup2,    $match_qual1, $match_qual2, $seqid2
);
my (
  @score1,      @score2, $cigar1, $cigar2, $ref_dna1, $ref_dna2, $first_mate,
  $second_mate, $seq_name
);
my @pairs = $sam->get_features_by_location(
  -type   => 'read_pair',
  -seq_id => "$chr",

  #  -start  => 500,
  #-end  => 80000000
);
my $num = 0;
my %snp_indel;
my %reads_type;
my %chain;
my %dobule_reads;
for my $pair (@pairs) {
  my %reads;
  $length = $pair->length;    # insert length
  ( $first_mate, $second_mate )
    = $pair->get_SeqFeatures;    #distinguish first and second reads
  $seqid1 = $first_mate->seq_id; #first reads chr_id
  if ( $first_mate && $second_mate ) {
    $start1 = $first_mate->start;     #first reads start mapping pos
    $start2 = $second_mate->start;    #second reads start mapping pos
    if ( ( $start1 > 0 ) && ( $start2 > 0 ) ) {
      $strand1 = $first_mate->strand;     #first reads forward or reverse
      $strand2 = $second_mate->strand;    #second reads forward or reverse
      $md1  = $first_mate->get_tag_values('MD');   #first reads mismatch base
      $md2  = $second_mate->get_tag_values('MD');  #second reads mismatch base
      $dup1 = $first_mate->get_tag_values('XA');   #first reads dup
      $dup2 = $second_mate->get_tag_values('XA');  #second reads dup
      $match_qual1 = $first_mate->qual;   #first reads mapping quality
      $match_qual2 = $second_mate->qual;  #second reads mapping quality
      @score1      = $first_mate->qscore; #first reads per-base quality scores
      @score2 = $second_mate->qscore;    #second reads per-base quality scores
      $cigar1 = $first_mate->cigar_str;  #first reads mapping
      $cigar2 = $second_mate->cigar_str; #second reads mapping
      $ref_dna1 = $first_mate->dna;
      $ref_dna2 = $second_mate->dna;
      $reads1   = $first_mate->query->dna;
      $reads2   = $second_mate->query->dna;
      $seq_name = $first_mate->name;
      my $md_cigar1;
      my $md_cigar2;

      if ( $md1
        && $cigar1
        && ( &filter_multialign_low_mq( $dup1, $match_qual1 ) ) ) {
        $md_cigar1 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md1, $cigar1, $reads1 ) )
          )
        );
      }
      else {
        next;
      }
      if ( $md2
        && $cigar2
        && ( &filter_multialign_low_mq( $dup2, $match_qual2 ) ) ) {
        $md_cigar2 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md2, $cigar2, $reads2 ) )
          )
        );
      }
      else {
        next;
      }
      my $index1;
      my $index2;
      my $reads1_length = length($reads1);
      my $reads2_length = length($reads2);
      if ( ( $reads1 =~ /^............ACT/ ) && ( $strand1 == 1 ) ) {
        $index1 = substr( $reads1, 0, 12 );
      }
      elsif ( ( $reads1 =~ /AGT............$/ ) && ( $strand1 == -1 ) ) {
        $index1 = reverse( substr( $reads1, $reads1_length - 12, 12 ) );
        $index1 =~ tr/ATCG/TAGC/;

        #print "reads1\t$strand1\t$index1\t$reads1\n";
      }
      else {
        print "reads1\t$strand1\t$reads1\n";
      }
      if ( ( $reads2 =~ /^............ACT/ ) && ( $strand2 == 1 ) ) {
        $index2 = substr( $reads2, 0, 12 );
      }
      elsif ( ( $reads2 =~ /AGT............$/ ) && ( $strand2 == -1 ) ) {
        $index2 = reverse( substr( $reads2, $reads2_length - 12, 12 ) );
        $index2 =~ tr/ATCG/TAGC/;

        #print "reads2\t$strand2\t$index2\t$reads2\n";
      }
      else {
        print "reads2\t$strand2\t$reads2\n";
      }
      if ( ($index1)
        && ($index2)
        && ( $index1 !~ /N/ )
        && ( $index2 !~ /N/ ) ) {
        my $start    = ( $start1 > $start2 ) ? $start2 : $start1;
        my $quality1 = ( sum @score1 ) / $reads1_length;
        my $quality2 = ( sum @score2 ) / $reads2_length;
        print $LOG1
          "$index1\t$index2\t$seqid1\t$start\t$length\t$quality1\t$quality2\t$start1\t$start2\t$md_cigar1\t$md_cigar2\t";
        print $LOG1
          "$seq_name\t$reads1_length\t$reads2_length\t$md1\t$md2\t$cigar1\t$cigar2\n";
      }
    }
  }
}

sub filter_no_mismatch {
  my ( $md, $cigar ) = @_;
  if ( ( $cigar !~ "D|I" ) && ( $md !~ "A|T|C|G" ) ) {
    print $LOG "filter_no_mismatch\t$cigar\t$md\n";
    return "0";
  }
  else {
    return "1";
  }
}

sub filter_multialign_low_mq {
  my ( $xa, $mq ) = @_;
  if ( $xa || ( $mq < 30 ) ) {
    print $LOG "filter_multialign_low_mq\t$mq\n";
    return "0";
  }
  else {
    return "1";
  }
}

sub merge_md_cigar {
  my ( $md, $cigar, $reads ) = @_;
  my @cigar1 = split /\d+/, $cigar;
  my @cigar2 = split /\D+/, $cigar;
  my @md1    = split /\d+/, $md;
  my @md2    = split /\D+/, $md;
  my $pos_reads = 0;
  my $sclip     = 0;
  my $ins_pos   = 0;
  my $ins_len   = 0;
  my @md_cigar;
  my $ins_alt;
  my $match_end = 0;

  for my $i ( 1 .. $#cigar1 ) {
    if ( ( $cigar1[$i] eq "S" ) && ( $i == 1 ) ) {
      push @md_cigar, $cigar2[0], "S";
      $pos_reads += $cigar2[0];
      $sclip = $cigar2[0];
    }
    elsif ( $cigar1[$i] eq "I" ) {
      $ins_alt = substr( $reads, $pos_reads, $cigar2[ $i - 1 ] );
      push @md_cigar, $pos_reads, "-/$ins_alt";
      $ins_pos = $pos_reads;
      $pos_reads += $cigar2[ $i - 1 ];
      $ins_len = $cigar2[ $i - 1 ];

      #    $match_end += $cigar2[ $i - 1 ];
    }
    elsif ( $cigar1[$i] eq "M" ) {
      $pos_reads += $cigar2[ $i - 1 ];
      $match_end += $cigar2[ $i - 1 ];
    }
  }
  $pos_reads = $sclip;
  my $indel_add = 0;
  for my $i ( 1 .. $#md1 ) {
    $pos_reads += $md2[ $i - 1 ];
    if ( ( $pos_reads >= $ins_pos )
      && ( $indel_add == 0 )
      && ( $ins_pos > 0 ) ) {
      $indel_add++;
      $pos_reads += $ins_len;
    }
    if ( $md1[$i] =~ /\^/ ) {
      $md1[$i] =~ s/\^//;
      push @md_cigar, $pos_reads, "$md1[$i]/-";
    }
    else {
      my $alt = substr( $reads, $pos_reads, 1 );
      push @md_cigar, $pos_reads, "$md1[$i]/$alt";
      $pos_reads++;
    }
  }
  my $last = max grep { $_ =~ /^\d+$/ } @md_cigar;
  push @md_cigar, $match_end;
  my $line = join "\,", @md_cigar;

  #my $read_len = length($reads);
  #print $LOG "$md\t$cigar\t$line\t$last\t$read_len\n";
  return $line;
}

sub sort_md_cigar {
  my ($line) = @_;
  my @new = split /\,/, $line;
  my $start = $new[0];
  my $tmp1;
  my $tmp2;
  for ( my $i = 2; $i <= $#new - 2; $i += 2 ) {
    if ( $new[$i] < $start ) {
      $tmp1          = $new[ $i - 2 ];
      $tmp2          = $new[ $i - 1 ];
      $new[ $i - 2 ] = $new[$i];
      $new[ $i - 1 ] = $new[ $i + 1 ];
      $new[$i]       = $tmp1;
      $new[ $i + 1 ] = $tmp2;
      $start         = $new[$i];
    }
    else {
      $new[$i]       = $new[$i];
      $new[ $i + 1 ] = $new[ $i + 1 ];
      $start         = $new[$i];
    }
  }
  my $line1 = join "\,", @new;
  return $line1;
}

sub filter_border {
  my ($md_cigar) = @_;
  my @md_cigar = split /\,/, $md_cigar;
  my @new;
  my $start = 0;
  my $end   = $md_cigar[-1];
  if ( ( $#md_cigar >= 1 ) && ( $md_cigar[1] eq "S" ) ) {
    $start += $md_cigar[0];
    $end   += $md_cigar[0];
  }
  for ( my $i = 0; $i <= $#md_cigar - 1; $i += 2 ) {
    if (
      !(
           ( $md_cigar[ $i + 1 ] )
        && ( $md_cigar[ $i + 1 ] ne "S" )
        && ( $md_cigar[$i] <= $start + 5 ) || ( $md_cigar[$i] >= $end - 5 )
      )
      ) {
      push @new, $md_cigar[$i], $md_cigar[ $i + 1 ];
    }
  }
  push @new, $md_cigar[-1];
  $md_cigar = join "\,", @new;
  return $md_cigar;
}

sub filter_many_mismatch {
  my ($md_cigar) = @_;
  my $count      = ( $md_cigar =~ tr/\/// );
  my $count_ind  = ( $md_cigar =~ tr/\-// );
  if ( $count_ind >= 2 ) {
    print $LOG "filter two or more indel $md_cigar\n";
    return "0";
  }
  elsif ( $count <= 3 ) {
    return $md_cigar;
  }
  else {
    print $LOG "filter more than three mismatch $md_cigar\n";
    return "0";
  }
}

sub prepare_quality_filter {
  my (@md_cigar) = @_;
  for ( my $i = 1; $i < $#md_cigar; $i += 2 ) {
    my $quality = 0;
    if ( ( $md_cigar[$i] ne "S" )
      && ( $md_cigar[$i] !~ "-" ) ) {    ##snp
      return (
        $md_cigar[ $i - 1 ] - 5,
        $md_cigar[ $i - 1 ],
        $md_cigar[ $i - 1 ] + 5,
        $md_cigar[$i]
      );
    }
    elsif ( $md_cigar[$i] =~ "\-\/" ) {    ##insertion
      return (
        $md_cigar[ $i - 1 ] - 5,
        $md_cigar[ $i - 1 ],
        $md_cigar[ $i - 1 ] + 5 + length( $md_cigar[$i] ) - 2,
        $md_cigar[$i]
      );
    }
    elsif ( $md_cigar[$i] =~ "\/\-" ) {    ##deletion
      return (
        $md_cigar[ $i - 1 ] - 5,
        $md_cigar[ $i - 1 ],
        $md_cigar[ $i - 1 ] + 4,
        $md_cigar[$i]
      );
    }
  }
}

sub pos_map {
  my ( $start1, $start2 ) = @_;
  my $reads_pos;
  if ( $start1 > $start2 ) {
    $reads_pos = "$start2.$start1";
  }
  else {
    $reads_pos = "$start1.$start2";
  }
  return "$reads_pos";
}

sub filter_low_quality {
  my (@quality) = @_;
  my $line = join "\t", @quality;
  my $num_q = 0;
  if ( $#quality == 10 ) {
    if ( $quality[5] > 20 ) {
      foreach (@quality) {
        if ( $_ > 20 ) {
          $num_q++;
        }
      }
      if ( $num_q >= 9 ) {
        return "1";
      }
      else {
        print $LOG "filter_low_quality:\n";
        return "0";
      }
    }
  }
  else {
    foreach (@quality) {
      if ( $_ > 20 ) {
        $num_q++;
      }
      elsif ( !$_ ) {
        print "err\t";
      }
    }
    if ( $num_q / @quality >= 0.8 ) {
      return "1";
    }
    else {
      print $LOG "filter_low_quality:\n";
      return "0";
    }
  }
}
