use strict;
open LIST, "$ARGV[0]";    ##cos的区段的bed文件
open IN,   "$ARGV[1]"
  ;    ##深度文件  第一列染色体  第二列位置  第三列深度
open OUT, ">$ARGV[1].bed";
open LOG, ">$ARGV[1].OUT";
my %hash;
my %chr;
my %dp;
my @sum;
my %hash1;
my $i = 0;
my $j;

#my $j=<LIST>;
#print OUT "$j";
while (<LIST>) {
  $i++;
  chomp;
  my @as = split;

  #		my @chr=split /chr/,$as[0];
  $hash1{ $as[0] } = 1;
  $chr{ $as[0] }{ $as[1] }{ $as[2] } = $_;
  if ( $sum[ $as[0] ] ) {
    $sum[ $as[0] ] = "$sum[$as[0]]\t$as[1]\t$as[2]";
  }
  else {
    $sum[ $as[0] ] = "$as[1]\t$as[2]";
  }
  $dp{ $as[0] }{ $as[1] }{ $as[2] }   = 0;
  $hash{ $as[0] }{ $as[1] }{ $as[2] } = 0;
}
my $sum  = 0;
my $sum1 = 0;
my $max;
my $min;
while (<IN>) {
  chomp;
  my @as = split;
  my @pos = split /\t/, $sum[ $as[0] ];
  for ( $i = 0; $i <= $#pos; $i += 2 ) {
    if ( ( $pos[$i] <= $as[1] ) && ( $pos[ $i + 1 ] >= $as[1] ) ) {
      $dp{ $as[0] }{ $pos[$i] }{ $pos[ $i + 1 ] } += $as[2];
      $hash{ $as[0] }{ $pos[$i] }{ $pos[ $i + 1 ] }++;
      $sum++;
      last;
    }
  }
}
print LOG "$sum\t$sum1\n";
foreach $i ( sort { $a <=> $b } ( keys %hash1 ) ) {
  my @pos = split /\t/, $sum[$i];
  for ( $j = 0; $j <= $#pos; $j += 2 ) {
    if ( ( $dp{$i}{ $pos[$j] }{ $pos[ $j + 1 ] } > 0 )
      && ( $hash{$i}{ $pos[$j] }{ $pos[ $j + 1 ] } > 0 ) ) {
      my $dp = int( $dp{$i}{ $pos[$j] }{ $pos[ $j + 1 ] }
          / $hash{$i}{ $pos[$j] }{ $pos[ $j + 1 ] } );
      print OUT "$chr{$i}{$pos[$j]}{$pos[$j+1]}\t$dp\n";
    }
    else {
      print OUT "$chr{$i}{$pos[$j]}{$pos[$j+1]}\t0\n";
    }
  }
}
