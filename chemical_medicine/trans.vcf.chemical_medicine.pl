use strict;
open VCF,                 "$ARGV[0]";
open TARGET_CHEMICAL,     "$ARGV[1]";
open TARGET_CHEMICAL_OUT, ">$ARGV[2].CHEMICAL.medicine.xls";

my %rs;
while (<VCF>) {
  chomp;
  my @line = split;
  my @frequency = split /\:/, $line[9];
  if    ( $frequency[2] >= 0.8 ) { $frequency[2] = 1; }
  elsif ( $frequency[2] >= 0.3 ) { $frequency[2] = 0.5; }
  else                           { $frequency[2] = 0; }
  $rs{ $line[0] }{ $line[1] }{ $line[3] }{ $line[4] } = $frequency[2];
}

my $name    = <TARGET_CHEMICAL>;
my @head    = split /\t/, $name;
my $raw_num = 6;
$name = join "\t", @head[ 0 .. $raw_num - 1 ];
print TARGET_CHEMICAL_OUT "$name\n";
my $geno;
while (<TARGET_CHEMICAL>) {
  my @line = split /\t/, $_;
  if (
    exists $rs{ $line[$raw_num] }{ $line[ $raw_num + 1 ] }
    { $line[ $raw_num + 2 ] }{ $line[ $raw_num + 3 ] } ) {
    if ( $rs{ $line[$raw_num] }{ $line[ $raw_num + 1 ] }
      { $line[ $raw_num + 2 ] }{ $line[ $raw_num + 3 ] }
      == $line[ $raw_num + 4 ] ) {
      $geno = join "\t", @line[ 0 .. $raw_num - 1 ];
      print TARGET_CHEMICAL_OUT "$geno\n";
    }
  }
  elsif ( $#line > $raw_num + 4 ) {
    my $true = 0;
    for ( my $i = $raw_num + 5; $i < $#line; $i = $i + 4 ) {
      if (
        exists $rs{ $line[$i] }{ $line[ $i + 1 ] }{ $line[ $i + 2 ] }
        { $line[ $i + 3 ] } ) {
        if ( $rs{ $line[$raw_num] }{ $line[ $raw_num + 1 ] }
          { $line[ $raw_num + 2 ] }{ $line[ $raw_num + 3 ] }
          == $line[ $raw_num + 4 ] ) {
          $geno = join "\t", @line[ 0 .. $raw_num - 1 ];
          print TARGET_CHEMICAL_OUT "$geno\n";
          $true++;
        }
      }
    }
    if ( $true == 0 ) {
      if ( $line[ ( $raw_num + 4 ) ] == 0 ) {
        $geno = join "\t", @line[ 0 .. $raw_num - 1 ];
        print TARGET_CHEMICAL_OUT "$geno\n";
      }
    }
  }
  else {
    if ( $line[ ( $raw_num + 4 ) ] == 0 ) {
      $geno = join "\t", @line[ 0 .. $raw_num - 1 ];
      print TARGET_CHEMICAL_OUT "$geno\n";
    }
  }
}
