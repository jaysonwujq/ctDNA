use strict;
use warnings;
use IO::All -utf8;

my $in         = shift;
my $target_bed = shift;

my @filter_header = qw(
left_chr
left_pos
left_strand
left_soft-clipped_reads
right_chr
right_pos
right_strand
right_soft-clipped_reads
SV_type
coverage_at_left_pos
coverage_at_right_pos
assembled_length_at_left_pos
assembled_length_at_right_pos
average_percent_identity_at_left_pos
percent_of_non-unique_mapping_reads_at_left_pos
average_percent_identity_at_right_pos
percent_of_non-unique_mapping_reads_at_right_pos
start_position_of_consensus_mapping_to_genome
starting_chromosome_of_consensus_mapping
position_of_the_genomic_mapping_of_consensus_starting_position
end_position_of_consensus_mapping_to_genome
ending_chromsome_of_consnesus_mapping
position_of_genomic_mapping_of_consensus_ending_posiiton
consensus_sequences
);
my $filter_header = join "\t", @filter_header;

sub sv_filter {
  my ( $percent_identity1, $percent_identity2, $depth1, $depth2 ) = @_;
  if ( ( $percent_identity1 >= 0.9 ) && ( $percent_identity2 >= 0.9 ) && ($depth1 >=4) && ($depth2>= 4) ) {
    return 1;
  }
  else {
    return 0;
  }
}

sub target_regio {
  my ( $chr, $pos ) = @_;
  my $target_all_pos = io($target_bed)->all;
  my @lines = split /\n/, $target_all_pos;
  foreach (@lines) {
    my @line = split /\t/, $_;
    if ( $chr > $line[0] ) {
      next;
    }
    elsif ( ( $chr eq $line[0] ) ) {
      if ( ( $line[1] - 50 <= $pos ) && ( $line[2] + 50 >= $pos ) ) {
        return 1;
      }
      else {
        next;
      }
    }
    elsif ( $chr < $line[0] ) {
      return 0;
    }
  }
  return 0;
}

my $sv = io($in)->all;
my @lines = split /\n/, $sv;
my @filter;
my @cnv;
push @filter, $filter_header;
foreach (@lines) {
  my @line = split /\t/, $_;
  if ( &sv_filter( $line[13], $line[15], $line[9], $line[10] ) ) {
    if ( ( &target_regio( $line[0], $line[1] ) )
      || ( &target_regio( $line[4], $line[5] ) ) ) {
      push @cnv,
        "$line[0]\t$line[1]\t$line[1]\tA\tG\thet\t.";
      push @cnv,
        "$line[4]\t$line[5]\t$line[5]\tA\tG\thet\t.";
      push @filter, $_;
    }
  }
}
my $cnv    = join "\n", @cnv;
my $filter = join "\n", @filter;
if (@filter > 1){
  $cnv > io("$in.filter.anno.txt");
  $filter > io("$in.filter");
}

