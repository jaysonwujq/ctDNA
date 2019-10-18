#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $chr );
use List::Util qw(max);
GetOptions(
  "in|i:s"  => \$in,
  "out|o:s" => \$out,

  #  "chr:s" => \$chr
);
use Bio::DB::Sam;
my $sam = Bio::DB::Sam->new(
  -bam   => "$in",
  -fasta => "/data1/software/b37/human_g1k_v37.fasta",
);
open OUT,  ">$out";
open LOG,  ">$out.filter";
open LOG1, ">$out.log";
my (
  $pair,    $length,      $seqid1,      $f_start, $start1,
  $start2,  $reads1,      $reads2,      $strand1, $strand2,
  $paired1, $paired2,     $md1,         $md2,     $dup1,
  $dup2,    $match_qual1, $match_qual2, $seqid2
);
my (
  @score1, @score2, $cigar1, $cigar2, $ref_dna1, $ref_dna2, $first_mate,
  $second_mate
);
my @pairs = $sam->get_features_by_location(
  -type => 'read_pair',

  #-seq_id => "$chr",
  #-start  => 500,
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
      my $md_cigar1;
      my $md_cigar2;

      if ( $md1
        && $cigar1
        && ( &filter_no_mismatch( $md1, $cigar1 ) )
        && ( &filter_multialign_low_mq( $dup1, $match_qual1 ) ) ) {
        $md_cigar1 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md1, $cigar1, $reads1 ) )
          )
        );
        print LOG1 "$md1\t$cigar1\t$md_cigar1\n";
        my @md_cigar1 = split /\,/, $md_cigar1;
        my $quality;
        my $start = 0;
        my $reads_pos;
        my @md_quality_pre = &prepare_quality_filter(@md_cigar1);
        for ( my $i = 0; $i < $#md_quality_pre; $i += 4 ) {
          $quality = 0;
          my @so
            = @score1[ $md_quality_pre[$i] .. $md_quality_pre[ $i + 2 ] ];

          #my $line = join "\t", @so;
          $quality = &filter_low_quality(@so);

          if ( $md_cigar1[1] eq "S" ) {
            $start = $md_cigar1[0];
          }
          if ( $quality == 1 ) {
            my $pos = $start1 + $md_quality_pre[ $i + 1 ] - $start;
            my $key = "$seqid1\t$pos\t$md_quality_pre[$i+3]";
            $snp_indel{$key}++;
            if ( !$start2 ) { $start2 = 0; }
            my $reads_pos = &pos_map( $start1, $start2 );
            if ( ( exists $reads_type{$key} )
              && ( $reads_type{$key} !~ /$reads_pos/ ) ) {
              $reads_type{$key} = "$reads_type{$key}\t$reads_pos";
              $chain{$key}{$strand1}++;
              $reads{$key}++;
            }
            elsif ( !exists $reads_type{$key} ) {
              $reads_type{$key} = "$reads_pos";
              $chain{$key}{$strand1}++;
              $reads{$key}++;
            }
          }
        }
      }

      if ( $md2
        && $cigar2
        && ( &filter_no_mismatch( $md2, $cigar2 ) )
        && ( &filter_multialign_low_mq( $dup2, $match_qual2 ) ) ) {
        $md_cigar2 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md2, $cigar2, $reads2 ) )
          )
        );
        print LOG1 "$md2\t$cigar2\t$md_cigar2\n";
        my @md_cigar2 = split /\,/, $md_cigar2;
        my $quality;
        my $start = 0;
        my $reads_pos;
        my @md_quality_pre = &prepare_quality_filter(@md_cigar2);
        for ( my $i = 0; $i < $#md_quality_pre; $i += 4 ) {
          $quality = 0;
          my @so
            = @score2[ $md_quality_pre[$i] .. $md_quality_pre[ $i + 2 ] ];
          $quality = &filter_low_quality(@so);
          if ( $md_cigar2[1] eq "S" ) {
            $start = $md_cigar2[0];
          }
          if ( $quality == 1 ) {
            my $pos = $start2 + $md_quality_pre[ $i + 1 ] - $start;
            my $key = "$seqid1\t$pos\t$md_quality_pre[$i+3]";
            $snp_indel{$key}++;
            if ( !$start1 ) { $start1 = 0; }
            my $reads_pos = &pos_map( $start1, $start2 );
            if ( ( exists $reads_type{$key} )
              && ( $reads_type{$key} !~ /$reads_pos/ ) ) {
              $reads_type{$key} = "$reads_type{$key}\t$reads_pos";
              $chain{$key}{$strand2}++;
            }
            elsif ( !exists $reads_type{$key} ) {
              $reads_type{$key} = "$reads_pos";
              $chain{$key}{$strand2}++;
            }
            elsif ( exists $reads{$key} ) {
              $chain{$key}{$strand2}++;
              $dobule_reads{$key}++;
            }
          }
        }
      }
    }
  }
}
foreach ( keys %snp_indel ) {
  my $type = split /\t/, $reads_type{$_};
  if ( !exists $chain{$_}{-1} )    { $chain{$_}{-1}    = 0; }
  if ( !exists $chain{$_}{1} )     { $chain{$_}{1}     = 0; }
  if ( !exists $dobule_reads{$_} ) { $dobule_reads{$_} = 0; }
  print OUT
    "$_\t$snp_indel{$_}\t$type\t$chain{$_}{-1}\t$chain{$_}{1}\t$dobule_reads{$_}\t$reads_type{$_}\n";
}

sub filter_no_mismatch {
  my ( $md, $cigar ) = @_;
  if ( ( $cigar !~ "D|I" ) && ( $md !~ "A|T|C|G" ) ) {
    print LOG "filter_no_mismatch\t$cigar\t$md\n";
    return "0";
  }
  else {
    return "1";
  }
}

sub filter_multialign_low_mq {
  my ( $xa, $mq ) = @_;
  if  ( $mq < 30 )  {
    print LOG "filter_multialign_low_mq\t$mq\n";
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
  #print LOG "$md\t$cigar\t$line\t$last\t$read_len\n";
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
  if ( $md_cigar[1] eq "S" ) {
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
    print LOG "filter two or more indel $md_cigar\n";
    return "0";
  }
  elsif ( $count <= 3 ) {
    return $md_cigar;
  }
  else {
    print LOG "filter more than three mismatch $md_cigar\n";
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
        print LOG "filter_low_quality:\n";
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
      print LOG "filter_low_quality:\n";
      return "0";
    }
  }
}
