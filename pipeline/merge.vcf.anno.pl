use strict;
open VCF,  "$ARGV[0]";
open ANNO, "$ARGV[1]";
open OUT,  ">$ARGV[1].merge";
my %hash;

while (<VCF>) {
  if ( $_ =~ "#" ) { next; }
  else {
    chomp;
    my @as1 = split;
    my @as2 = split /\:/, $as1[9];
    my @as3 = split /\,/, $as2[1];

    #		  $hash{$as1[0]}{$as1[1]}{$as1[3]}{$as1[4]}=""
    my $anno = <ANNO>;
    chomp $anno;
    my @as4 = split /\t/, $anno;
    my $line1 = join "\t", @as4[ 0 .. 4 ];
    my $line2 = join "\t", @as4[ 22 .. $#as4 - 6 ];

    #if ($as4[-1]!~"COS"){$as4[-1]="NA";}
    print OUT
      "$line1\t$as1[2]\t$as1[-1]\t$as1[5]\t$as1[-3]\t$as1[-2]\t$as4[8]\t$as4[7]\t$as4[10]\t$as4[11]\t$as4[21]\t$as3[0]\t$as3[1]\t$as2[2]\t$as2[3]\t$as2[4]\t$line2\n";
  }
}
