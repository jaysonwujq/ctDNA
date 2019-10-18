use  warnings;
my $sample = $ARGV[0];

my $in = $ARGV[1];
my %hash;
my @line;
my $header = "sample_ID\ttotal_mapped_reads\ttarget_reads\ttarget_ratio";
chomp $sample;
open IN, "$in/$sample.sort.bam.flagstat";
<IN>;
<IN>;
my $line3 = <IN>;
<IN>;
my $line5 = <IN>;
my $line;
if   ( $line3 =~ /map/ ) { $line = $line3; }
else                     { $line = $line5; }
my @as = split /\s+/, $line;
$line[0] = $as[0];
close IN;
open IN, "$in/$sample.target.bam.flagstat";
<IN>;
<IN>;
$line3 = <IN>;
<IN>;
$line5 = <IN>;
if   ( $line3 =~ /map/ ) { $line = $line3; }
else                     { $line = $line5; }
@as = split /\s+/, $line;
$line[1] = $as[0];
close IN;
$line[2] = $line[1] / $line[0];
my $lsi = join "\t", @line;
my $email = 'liubei@cheerlandgroup.com,zhaomingyu@cheerlandgroup.com';
my $message = "total_mapped_reads:\t$line[0]\ntarget_reads:\t$line[1]\ntarget_ratio:\t$line[2]\n";
`echo "$sample\n$message"| /usr/bin/mail -s "$sample" $email`;
