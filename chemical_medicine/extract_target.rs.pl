open TARGET,
  "/data1/workdir/liubei/program/ctDNA/data_base/target_bed/cosmic.list.bed.old2";
open RS,         "$ARGV[0]";
open TARGET_ALF, ">$ARGV[0].target";

my %target;
my @pos;
while (<TARGET>) {
  chomp;
  @pos = split;
  foreach ( $pos[1] - 100 .. $pos[2] + 100 ) {
    $target{ $pos[0] }{$_} = 1;
  }
}
my $name = <RS>;
chomp $name;
my @tmp = split /\t/, $name;
$name = join "\t", @tmp[ 0 .. $#tmp - 4 ];
print TARGET_ALF "$name\tchr\tpos\tref\talt\tgeno\n";

my @rs;
while (<RS>) {
  chomp;
  @rs = split /\t/, $_;
  if ( exists $target{ $rs[-4] }{ $rs[-3] } ) {
    if ( $rs[-1] =~ "," ) {
      my @alt = split /\,/, $rs[-1];
      foreach my $k ( 0 .. $#alt ) {
        if ( $rs[3] ne "$rs[-2]$rs[-2]" ) {
          if ( $rs[3] =~ "$alt[$k]" ) {
            $rs[-1] = $alt[$k];
          }
        }
      }
    }
    if    ( $rs[3] eq "$rs[-2]$rs[-2]" ) { $rs[ $#rs + 1 ] = 0; }
    elsif ( $rs[3] eq "$rs[-1]$rs[-1]" ) { $rs[ $#rs + 1 ] = 1; }
    elsif ( ( $rs[3] eq "$rs[-1]$rs[-2]" ) || ( $rs[3] eq "$rs[-2]$rs[-1]" ) )
    {
      $rs[ $#rs + 1 ] = 0.5;
    }
    elsif ( $rs[3] eq "del/del" ) { $rs[ $#rs + 1 ] = 1; }
    elsif ( $rs[3] =~ "del" )     { $rs[ $#rs + 1 ] = 0.5; }
    elsif ( $rs[3] !~ "6|7" )     { $rs[ $#rs + 1 ] = 0; }
    my $all = join "\t", @rs;
    print TARGET_ALF "$all\n";
  }
  elsif ( $rs[3] =~ "6|7" ) {
    if ( ( $rs[3] =~ "6" ) && ( $rs[3] =~ "7" ) ) { $rs[ $#rs + 1 ] = 1; }
    elsif ( $rs[3] =~ "7" ) { $rs[ $#rs + 1 ] = 0; }
    elsif ( $rs[3] =~ "6" ) { $rs[ $#rs + 1 ] = 0.5; }
    @rs[ -5 .. -2 ] = ( "2", "234668880", "ATA", "A" );
    my $all = join "\t", @rs;
    print TARGET_ALF "$all";
    for ( my $i = 234668881; $i <= 234668891; $i = $i + 2 ) {
      print TARGET_ALF "\t2\t$i\tTAT\tT";
    }
    for ( my $i = 234668882; $i <= 234668892; $i = $i + 2 ) {
      print TARGET_ALF "\t2\t$i\tATA\tA";
    }
    print TARGET_ALF "\t2\t234668893\tTAA\tA\n";
  }
}
