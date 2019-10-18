use strict;
use warnings;
my $sv  =shift;
my $cds = shift;
open my $SV, "<", "$sv";
open my $OUT, ">", "$sv.exon.intron";
open my $OUT1,">", "$sv.filtered";
my $header = <$SV>;
chomp $header;
my @headers = split /\t/, $header;
$headers[2] = "$headers[2]\tgene_strand\tintro_or_exon_num";
$headers[8] = "$headers[8]\tgene_strand\tintro_or_exon_num";
$header = join "\t", @headers;
print $OUT  "$header\n";
print $OUT1 "type\t$header\n";

while (<$SV>){
  chomp;
  my @line = split /\t/,$_;
  my $tmp1 = `grep -P \"\\t$line[3]\$\" $cds`;
  my $tmp2 = `grep -P \"\\t$line[9]\$\" $cds`;
  my $gene1 = &num_exon_intro($line[2], $line[1], $tmp1);
  my $gene2 = &num_exon_intro($line[8], $line[7], $tmp2);
  $line[2] = "$line[2]\t$gene1";
  $line[8] = "$line[8]\t$gene2";
  my $result =join "\t", @line;
  print $OUT "$result\n";
  if ( ($gene1 eq $gene2) && ($line[2] =~ /intronic/) &&($line[8] =~ /intronic/) && ($line[3] eq $line[9]) ){
    print $OUT1 "filtered-intronic-same\t$result\n";
  }
  elsif (($line[13] <=  300) && ($line[14] <= 300)){
    print $OUT1 "filtered-low_coverage\t$result\n";
  }
  elsif (($line[2] =~ /intergenic/ ) ||($line[8] =~ /intergenic/)){
    print $OUT1 "filtered-intergenic\t$result\n";
  }
  else{
    my $strand1 = &strand($line[4]);
    my $strand2 = &strand($tmp1);
    my $strand3 = &strand($line[10]);
    my $strand4 = &strand($tmp2);
    if ( ($strand1*$strand3) == ($strand2*$strand4) ){
      print $OUT1 "ok\t$result\n";
    }
    else{
      print $OUT1 "filtered-strand\t$result\n";
    }
  }
}

sub strand{
  my $strand = $_[0];
  if ($strand =~ /\+/){
    return 1;
  }
  else{
    return -1;
  }
}

sub num_exon_intro{
  ## get the postion of the gene 
  my ($type, $pos, $gene) = @_;
  my @lines = split /\n/, $gene;
  if ($type eq "exonic"){
    foreach (@lines){
      my @line = split  /\t/, $_;
      if ( ($line[1] < $pos) && ($line[2] > $pos) ){
        return "$line[3]\t$line[6]";
      }
    }
  }
  elsif (($type eq "intronic") && ($gene =~ /\+/)){
    foreach (@lines){
      my @line = split  /\t/, $_;
      if ($pos < $line[1]){
        my $num = $line[6] - 1;
        return "$line[3]\t$num";
      }
    }
  }
  elsif (($type eq "intronic") && ($gene !~ /\+/)){
    foreach (@lines){
      my @line = split  /\t/, $_;
      if ($pos > $line[1]){
        my $num = $line[6] - 1;
        return "$line[3]\t$num";
      }
    }
  }else{
    return "NA\tNA";
  }
}

