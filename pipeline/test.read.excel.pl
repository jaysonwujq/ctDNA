#!/bin/usr/perl -w

use Text::Iconv;
use strict;
use Getopt::Long;
my ( $in, $type, $out );
GetOptions(
  "in|i:s"  => \$in,
  "type:s"  => \$type,
  "out|o:s" => \$out
);
my $converter = Text::Iconv->new( "utf-8", "windows-1251" );
open OUT, ">$out.$type.split";

# Text::Iconv is not really required.
#  # This can be any object with the convert method. Or nothing.
#
use Spreadsheet::XLSX;
my @all;
my $excel = Spreadsheet::XLSX->new( "$in", $converter );
foreach my $sheet ( @{ $excel->{Worksheet} } ) {
  printf( "Sheet: %s\n", $sheet->{Name} );
  if ( $sheet->{Name} =~ "sample|Sample" ) {
    $sheet->{MaxRow} ||= $sheet->{MinRow};
    foreach my $row ( $sheet->{MinRow} .. $sheet->{MaxRow} ) {
      my @result;
      $sheet->{MaxCol} ||= $sheet->{MinCol};
      foreach my $col ( $sheet->{MinCol} .. $sheet->{MaxCol} ) {
        my $cell = $sheet->{Cells}[$row][$col];
        if ($cell) {
          $result[$col] = $cell->{Val};
        }
      }
      my $row1 = join "\t", @result;
      $all[$row] = $row1;
    }
  }
}
my @name;
my @use = ( "sample", "project_type", "index1_seq", "index2_seq" );
my @index;
foreach (@all) {
  if ( $_ =~ "sample" ) {
    @name = split /\t/, $_;
    foreach my $i ( 0 .. $#name ) {
      foreach my $j ( 0 .. $#use ) {
        if ( $name[$i] eq "$use[$j]" ) {
          $index[$j] = $i;
        }
      }
    }
  }
  elsif ( $_ =~ /\t/ ) {
    my @line = split /\t/, $_;
    my @out;
    print "$_\n";
    if ( $line[ $index[1] ] eq "$type" ) {
      foreach my $i ( 0 .. $#line ) {
        foreach my $j ( 0 .. $#index ) {
          if ( ( $i == $index[$j] ) && ( $line[$i] =~ /\S+/ ) ) {
            if ( ( $line[$i] ne "index2_seq" )
              && ( $line[$i] ne "$type" ) ) {
              $line[$i] =~ s/\%|\s//;
              push @out, $line[$i];
            }
          }
        }
      }

      my $line = join "\t", @out;
      print OUT "$line\n";
    }
  }
}
