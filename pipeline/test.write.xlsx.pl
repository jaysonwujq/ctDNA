use strict;
use Excel::Writer::XLSX;
open IN, '<:utf8', $ARGV[0];
my $out = $ARGV[1];

# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new("$out");

# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Add and define a format
#$format = $workbook->add_format();
#$format->set_bold();
#$format->set_color( 'red' );
#$format->set_align( 'center' );
# Write a formatted and unformatted string, row and column notation.
#my $col = $row = 0;
#$worksheet->write( $row, $col, 'Hi Excel!', $format );
#$worksheet->write( 1, $col, 'Hi Excel!' );
# Write a number and a formula using A1 notation
#$worksheet->write( 'A3', 1.2345 );
#$worksheet->write( 'A4', '=SIN(PI()/4)' );
my $i = 0;
while (<IN>) {
  chomp;
  $i++;
  my @line = split /\t/, $_;

  #  print "$_\n";
  $worksheet->write( "A$i", \@line );
}
