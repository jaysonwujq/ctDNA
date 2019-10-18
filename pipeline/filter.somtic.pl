open VCF,  "$ARGV[0]";
open ANNO, "$ARGV[1]";
open OUT,  ">$ARGV[2]";

while (<VCF>) {
  if ( $_ =~ "#" ) { next; }
  else {
    chomp;
    my @vcf = split /\t/, $_;
    my $anno = <ANNO>;
    chomp $anno;
    my @anno = split /\t/, $anno;
    if ( ( ( $anno[15] + $anno[16] ) >= 30 )
      && ( $anno[16] >= 5 )
      && ( $anno[17] >= 0.002 ) ) {
      print OUT "$anno\n";
    }
  }
}
