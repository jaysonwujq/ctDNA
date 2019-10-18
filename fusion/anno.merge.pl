use strict;
use warnings;

open my $SV, "<", "$ARGV[0]";
open my $ANNO, "<", "$ARGV[1]";
open my $OUT, ">", "$ARGV[0].ANNO";

my $header = <$SV>;
chomp $header;
my @headers = split /\t/, $header;
$headers[1] = "$headers[1]\tpos_of_gene\tgene";
$headers[5] = "$headers[5]\tpos_of_gene\tgene";
$header = join "\t", @headers;
print $OUT "$header\n";

while (<$SV>){
  chomp;
  my @line = split /\t/,$_;
  my $line1 = <$ANNO>;
  my $line2 = <$ANNO>;
  my $anno1 = (split /\t/,$line1)[0];
  my $anno2 = (split /\t/,$line1)[1];
  my $anno3 = (split /\t/,$line2)[0];
  my $anno4 = (split /\t/,$line2)[1];
  $line[1] = "$line[1]\t$anno1\t$anno2";
  $line[5] = "$line[5]\t$anno3\t$anno4";
  my $result  = join "\t",@line;
  print $OUT "$result\n";
}

