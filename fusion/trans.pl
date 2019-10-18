use strict;

open my $BAM, "samtools view $ARGV[0]|";
my $out  = $ARGV[1];
my $name = $ARGV[2];

my $alk  = "2 29415640 30144432";
my $ros1 = "6 117609463 117747018";
my $ret  = "10 43572475 43625799";

open my $ALK,  ">", "$out/$name.ALK";
open my $ROS1, ">", "$out/$name.ROS1";
open my $RET,  ">", "$out/$name.RET";

my %reads;
my %alk;
my %ros1;
my %ret;
my %tmp;

my @alk = split /\s+/, $alk;
foreach ( $alk[1] .. $alk[2] ) {
  $tmp{ $alk[0] }{$_} = "ALK";
}

my @ros1 = split /\s+/, $ros1;
foreach ( $ros1[1] .. $ros1[2] ) {
  $tmp{ $ros1[0] }{$_} = "ROS1";
}

my @ret = split /\s+/, $ret;
foreach ( $ret[1] .. $ret[2] ) {
  $tmp{ $ret[0] }{$_} = "RET";
}

while (<$BAM>) {
  chomp;
  my @lines = split /\t/, $_;
  if ( exists $tmp{ $lines[2] }{ $lines[3] } ) {
    $reads{ $lines[0] } = $tmp{ $lines[2] }{ $lines[3] };
  }
}

open my $BAM, "samtools view $ARGV[0]|";
while (<$BAM>) {
  chomp;
  my @lines = split /\t/, $_;
  if ( exists $reads{ $lines[0] } ) {
    if ( $reads{ $lines[0] } eq "ALK" ) {
      say $ALK "$_";
    }
    elsif ( $reads{ $lines[0] } eq "RET" ) {
      say $RET $_;
    }
    elsif ( $reads{ $lines[0] } eq "ROS1" ) {
      say $ROS1 $_;
    }
  }
}
