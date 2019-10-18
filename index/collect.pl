use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $out, $list, $sam );
GetOptions(
  "out|o:s"  => \$out,
  "list|l:s" => \$list,
  "sam:s"    => \$sam,
);
open my $OUT1,  ">", "$out.1";
open my $OUT2,  ">", "$out.2";
open my $OUT3,  ">", "$out.3";
open my $OUT4,  ">", "$out.4";

open my $BAM,  "<", "$sam";
open my $LIST, "<", "$list";
my %reads;
while (<$LIST>) {
  chomp;
  my @tmp = split /\s+/, $_;
  my $reads_name = $tmp[0];
  my $reads_num  = $tmp[1];
  $reads{$reads_name} = $reads_num ;
}

while (<$BAM>) {
  chomp;
  if ( $_ =~ /^@/ ) {
    print $OUT1 "$_\n";
    print $OUT2 "$_\n";
    print $OUT3 "$_\n";
    print $OUT4 "$_\n";
  }
  else {
    my $seq_id = ( split /\t/, $_ )[0];
    if ( exists $reads{$seq_id} ) {
      if ( $reads{$seq_id} >= 1){
        print $OUT1 "$_\n";
      }
      if ( $reads{$seq_id} >= 2){
        print $OUT2 "$_\n";
      }
      if ( $reads{$seq_id} >= 3){
        print $OUT3 "$_\n";
      }
      if ( $reads{$seq_id} >= 4){
        print $OUT4 "$_\n";
      }
    }
  }
}
