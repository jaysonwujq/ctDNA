open LIST, "$ARGV[0]";
my $outdir = $ARGV[1];

#my $indir=$ARGV[2];
`mkdir -p $outdir`;
open TARGET,
  "/data1/workdir/liubei/progrem/cosmic/pipline/cosmic.list.bed.old2";
use Statistics::Basic qw(:all nofill);
my %hash;
my %median;
my %target;
my $pos_all = 0;
open LOG, ">$outdir/target.pos";

while (<TARGET>) {
  chomp;
  my @pos = split;
  foreach my $i ( $pos[1] .. $pos[2] ) {
    $pos_all++;
    $target{ $pos[0] }{$i} = 1;
    print LOG "$pos_all\t$pos[0]\t$i\n";
  }
}
close TARGET;
while (<LIST>) {
  chomp;
  my $sample = $_;
  open IN, "$outdir/$sample.target.bam.depth";
  my $line = 0;
  my @dp;
  while (<IN>) {
    chomp;
    my @as = split;
    if ( exists $target{ $as[0] }{ $as[1] } ) {
      $hash{ $as[0] }{ $as[1] }{$sample} = "$as[2]";
      push @dp, $as[2];
    }
  }
  $median{$sample} = median(@dp);
  print "$median{$sample}\t";
}
my @smaples;
foreach my $sample ( sort keys %median ) {
  push @samples, $sample;
}
my $samples = join "\t", @samples;
my $ratio;

open TARGET,
  "/data1/workdir/liubei/program/ctDNA/cnv/cosmic.list.bed.old2.OUT.sort2";
my $num = 0;

#open OUT,">all.dp";
print OUT "$samples\n";
my $before = 0;
while (<TARGET>) {
  chomp;
  $num++;
  my @pos = split;
  if ( $pos[4] ne $before ) {
    open OUT, ">$outdir/$pos[4]";
    $before = $pos[4];
    print OUT "$samples\n";
  }

  #  open OUT,">$outdir/section$num";
  #  print OUT "$samples\n";
  foreach my $i ( $pos[1] .. $pos[2] ) {
    my @ratio;
    foreach my $sample ( sort keys %median ) {
      if ( exists $hash{ $pos[0] }{$i}{$sample} ) {
        $ratio = $hash{ $pos[0] }{$i}{$sample} / $median{$sample} * 100;
      }
      else {
        $ratio = "NA";
      }
      push @ratio, $ratio;
    }
    my $ratios = join "\t", @ratio;
    print OUT "$ratios\n";
  }
}
