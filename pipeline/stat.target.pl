open LIST, "$ARGV[0]";
my $in = $ARGV[1];
open OUT, ">$in/target.ratio";
my %hash;
my @line;
print OUT "sample_ID\ttotal_mapped_reads\ttarget_reads\ttarget_ratio\n";
while (<LIST>) {
  chomp;
  open IN, "$in/$_.sort.bam.flagstat";
  <IN>;
  <IN>;
  my $line3 = <IN>;
  <IN>;
  my $line5 = <IN>;
  my $line;
  if   ( $line3 =~ /map/ ) { $line = $line3; }
  else                     { $line = $line5; }
  my @as = split /\s+/, $line;

  #		  print "$_\t$as[0]\t";
  $line[0] = $as[0];
  close IN;
  open IN, "$in/$_.target.bam.flagstat";

  #print "$_.target.bam.flagstat\t";
  <IN>;
  <IN>;
  my $line3 = <IN>;
  <IN>;
  my $line;
  my $line5 = <IN>;
  if   ( $line3 =~ /map/ ) { $line = $line3; }
  else                     { $line = $line5; }
  my @as = split /\s+/, $line;

  #		  print "$as[0]\t";
  $line[1] = $as[0];
  close IN;
  $line[2] = $line[1] / $line[0];
  my $lsi = join "\t", @line;
  print OUT "$_\t$lsi\n";
}

