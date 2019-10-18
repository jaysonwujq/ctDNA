#!/bin/usr/perl -w

use Text::Iconv;
use strict;
use Getopt::Long;
my ( $in, $type, $out );
GetOptions(
  "in|i:s" => \$in,
);

#my $converter = Text::Iconv->new( "windows-1251","windows-1251");
my $converter = Text::Iconv->new( "utf8", "utf8" );

# Text::Iconv is not really required.
#  # This can be any object with the convert method. Or nothing.
#
use Spreadsheet::XLSX;
my @all;
my $excel = Spreadsheet::XLSX->new( "$in", $converter );
foreach my $sheet ( @{ $excel->{Worksheet} } ) {
  printf( "Sheet: %s\n", $sheet->{Name} );
  my $out1 = $sheet->{Name};
  open my $OUT1, ">", "$out1";
  $sheet->{MaxRow} ||= $sheet->{MinRow};
  foreach my $row ( $sheet->{MinRow} .. $sheet->{MaxRow} ) {
    my @result;
    $sheet->{MaxCol} ||= $sheet->{MinCol};
    foreach my $col ( $sheet->{MinCol} .. $sheet->{MaxCol} ) {
      my $cell = $sheet->{Cells}[$row][$col];
      print $OUT1 "$cell->{Val}\t";
    }
    print $OUT1 "\n";
  }
}
