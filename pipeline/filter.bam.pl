#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $chr, $in1 );
use List::Util qw(max);
GetOptions(
  "in|i:s"  => \$in,
  "bam:s"   => \$in1,
  "out|o:s" => \$out,

  #  "chr:s" => \$chr
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
  $pair,    $length,      $seqid1,      $f_start, $start1,
  $start2,  $reads1,      $reads2,      $strand1, $strand2,
  $paired1, $paired2,     $md1,         $md2,     $dup1,
  $dup2,    $match_qual1, $match_qual2, $seqid2,  $seq_name,
  $nm1,     $nm2
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
my %count;
my %read1;
my %read2;
for my $pair (@pairs) {
  my %reads;
  $length = $pair->length;    # insert length
  ( $first_mate, $second_mate )
    = $pair->get_SeqFeatures;    #distinguish first and second reads

  $seqid1 = $first_mate->seq_id; #first reads chr_id
  if ( $first_mate && $second_mate ) {
    $start1 = $first_mate->start;     #first reads start mapping pos
    $start2 = $second_mate->start;    #second reads start mapping pos
    if ( ($start1) && ($start2) && ( $start1 > 0 ) && ( $start2 > 0 ) ) {
      $strand1 = $first_mate->strand;     #first reads forward or reverse
      $strand2 = $second_mate->strand;    #second reads forward or reverse
      $md1  = $first_mate->get_tag_values('MD');   #first reads mismatch base
      $md2  = $second_mate->get_tag_values('MD');  #second reads mismatch base
      $dup1 = $first_mate->get_tag_values('XA');   #first reads dup
      $dup2 = $second_mate->get_tag_values('XA');  #second reads dup
      $nm1  = $first_mate->get_tag_values('NM');   #First reads mismatch base numbers
      $nm2  = $first_mate->get_tag_values('NM');   #First reads mismatch base numbers
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
      $seq_name = $first_mate->name;
      if ( $md1
        && $cigar1
        && ( &slip_map( $md1, $cigar1 ) )
        && ( &filter_multialign_low_mq(  $match_qual1 ) ) ) {
        $md_cigar1 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md1, $cigar1, $reads1 ) )
          )
        );
        if ($md_cigar1) {
          $read1{$seq_name} = 1;
        }
        print $LOG1 "$seqid1\t$start1\t$start2\t$md_cigar1\n";
      }

      if ( $md2
        && $cigar2
        && ( &slip_map( $md2, $cigar2 ) )
        && ( &filter_multialign_low_mq(  $match_qual2 ) ) ) {
        $md_cigar2 = &filter_many_mismatch(
          &sort_md_cigar(
            &filter_border( &merge_md_cigar( $md2, $cigar2, $reads2 ) )
          )
        );
        if ($md_cigar2) {
          $read2{$seq_name} = 1;
        }
        print $LOG1 "$seqid1\t$start1\t$start2\n";
        my @md_cigar2 = split /\,/, $md_cigar2;
        print $LOG1 "$seqid1\t$start1\t$start2\t$md_cigar2\n";
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
    if ($tmp[1] < 128){
      if (exists $read1{$key}){
        print $OUT "$_\n";
      }
    }
    else{
      if (exists $read2{$key}){
        print $OUT "$_\n";
      }
    }
  }
}

sub slip_map {
  my ( $md, $cigar ) = @_;
  return "1";
}

sub filter_multialign_low_mq {
  my ( $mq ) = $_[0];
  if (  $mq < 30  ) {
    print $LOG "filter_low_mq\t$mq\n";
    return "0";
  }
  else {
    return "1";
  }
}
sub filter_multialign {
  my ( $xa1, $xa2, $nm1, $nm2 ) = @_;
  if ( $xa1 && $xa2 ) {
    my %class;
    my $xa1_nm1 = (split /\,|\;/, $xa1)[3];
    my $xa2_nm2 = (split /\,|\;/, $xa2)[3];
    my @nms     = ($nm1, $nm2, $xa1_nm1, $xa2_nm2);
    my @names   = qw (nm1, nm2, xa1_nm1, xa2_nm2);
    foreach my $i (0..$#nms){
      if ($ == 0 ){
        $class{$names[$i]}=1;
      }elsif ($_ <= 2 ){
        $class{$names[$i]}=2;
      }else{
        $class{$names[$i]}=3;
      }
    }
    if( ($class{nm1} <= $class{xa1_nm1}) && ($class{nm2} <= $class{xa2_nm2}) ){
      print $LOG "filter_multialign\t$xa1, $xa2, $nm1, $nm2\n";
      return "0";
    }
    else{
      print $LOG "$xa1, $xa2, $nm1, $nm2\n";
      return "1";
    }
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
  if ( $#md_cigar == 0 ) { return $md_cigar; }
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
  if ( $count_ind > 2 ) {
    print $LOG "filter_three_or_more_indel $md_cigar\n";
    return "0";
  }
  elsif ( $count <= 5 ) {
    return $md_cigar;
  }
  else {
    print $LOG "filter_more_than_five_mismatch $md_cigar\n";
    return "0";
  }
}

