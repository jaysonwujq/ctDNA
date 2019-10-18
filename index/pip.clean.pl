#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $list, $target, $var, $type, $fa, $anno );
GetOptions(
  "in:s"     => \$in,
  "out:s"    => \$out,
  "target:s" => \$target,
  "type:s"   => \$type,
  "fa:s"     => \$fa,
  "anno:s"   => \$anno,
);
`mkdir -p $out $out/shell $out/clean $out/bam $out/vcf $out/stat $out/bam.stat $out/anno $out/amplication $out/fusion $out/result $out/DP`;
`mkdir -p $out/release $out/release/bam $out/release/result`;
my $bin      = $Bin;
my $samtools = "/data1/workdir/liubei/program/WES/samtools";
if ( !$fa )   { $fa   = "/data1/software/b37/human_g1k_v37.fasta"; }
if ( !$anno ) { $anno = "/data1/software/annovar/"; }
open my $SNP_SH, ">", "$out/shell/raw_to_vcf.sh";
open my $DO_SH,  ">", "$out/shell/do.qsub.sh";

#print SPLIT_SH "perl $bin/test.read.excel.pl -in $list -type ctDNA -out $out/list";
#`perl $bin/test.read.excel.pl -in $list -type ctDNA -out $out/list`;
open my $LIST, "<", "$in/stats/list";
open my $TYPE, "<", "$type";
my $i = 0;
my $tumor_type;
my %case;
my %control;
my %sample;

while (<$LIST>) {
  chomp;
  open my $SH1,       ">", "$out/shell/$_.sh";
  open my $SHBWA,     ">", "$out/shell/$_.bwa.sh";
  open my $DOBWA,     ">", "$out/shell/do.$_.bwa.sh";
  open my $SHBWASTAT, ">", "$out/shell/$_.bwa.stat.sh";
  open my $SNP,       ">", "$out/shell/$_.snpcalling.sh";
  open my $ANNO,      ">", "$out/shell/$_.anno.sh";
  $i++;
  my $type_real = <$TYPE>;
  my @type = split /\t/, $type_real;
  if    ( $type[0] == 1 ) { $tumor_type = "all.RS.target"; }
  elsif ( $type[0] == 2 ) { $tumor_type = "lung.RS.target"; }
  elsif ( $type[0] == 3 ) { $tumor_type = "jzc.RS.target"; }

  if ( $type[2] == 0 ) {
    $control{ $type[1] } = 1;
    $sample{ $type[1] }{control} = $_;
  }
  else { $case{ $type[1] } = 1; $sample{ $type[1] }{case} = $_ }

  #my $clean=&raw_to_clean($_);
  #  print SNP_SH "$clean";
  my $bwa     = &mapping($_);
  my $rmdup   = &rmdup($_);
  my $realign = &realign($_);
  print $SHBWA "$bwa\n$rmdup\n$realign\n";
  my $bwa_stat    = &bam_stat($_);
  my $snp_calling = &snp_calling( $_, $tumor_type );
  my $anno_result = &anno( $_, $tumor_type );
  print $SHBWASTAT "$bwa_stat\n";
  print $SNP "$snp_calling\nsh $out/shell/$_.bwa.stat.sh\n";
  print $ANNO "$anno_result\n";

  #  print $SNP_SH "sh $out/shell/$_.sh\n";
  print $DOBWA "sh $out/shell/$_.bwa.sh\n";

  #  print $SH1 "sh $out/shell/$_.bwa.stat.sh\n";
  #print $SH1 "sh $out/shell/$_.snpcalling.sh\n";
  #print $SH1 "sh $out/shell/$_.anno.sh\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=8,mem=30G,walltime=1000:00:00 --maxjob 50 $out/shell/do.$_.bwa.sh --convert n --jobprefix ct.$i.bwa --queue batch \n";

  #print $SH1 "$bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=1,mem=10G,walltime=1000:00:00 --maxjob 50 $out/shell/$_.bwa.stat.sh\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 1  --resource nodes=1:ppn=3,mem=15G,walltime=1000:00:00 --maxjob 50 $out/shell/$_.snpcalling.sh --convert n --jobprefix ct.$i.snp --queue all\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=1,mem=10G,walltime=1000:00:00 --maxjob 50 $out/shell/$_.anno.sh --convert n  --jobprefix ct.$i.anno --queue batch\n";
  close $SH1;
  print $SNP_SH "sh $out/shell/$_.sh&\n";
}

#print $DO_SH "perl $bin/qsub-pbs.pl --lines 1  --resource nodes=1:ppn=4,mem=15G,walltime=1000:00:00 --maxjob 50 $out/shell/raw_to_vcf.sh \n";
# print $DO_SH "perl $bin/stat.target.pl $out/list $out/bam.stat\n";
#print $DO_SH "sh $out/shell/raw_to_vcf.sh\n";
print $DO_SH
  "perl $bin/target.dp.pl $in/stats/list $target $out/DP/all.dp.stat $out/DP\n";
print $DO_SH "perl $bin/stat.target.pl $in/stats/list $out/bam.stat\n";
print $DO_SH
  "perl $bin/merge.all.stat.pl $in/stats/clean.stat $out/bam.stat/target.ratio $out/DP/all.dp.stat $out/all.stat\n";
print $DO_SH
  "perl $bin/test.write.xlsx.pl $out/all.stat $out/all.stat.xlsx\n";
print $DO_SH "ln -sf $out/all.stat.xlsx $out/release/result/.\n";

#print $DO_SH "perl /data1/workdir/liubei/program/cosmic/pipline/new/qsub-pbs.pl --queue batch --lines 1  --resource nodes=1:ppn=4,mem=25G,walltime=1000:00:00   --maxjob 50 $out/shell/raw_to_vcf.sh\n";
#print DO_SH "cd $out/release\n python -m SimpleHTTPServer 5004 &\n";
#print $DO_SH "echo \"Hello, $out ctDNA data snp and indel calling  was done ; You can go to http://192.168.1.205:5004/ or $out/release  for the bam and result excel\"| mail -s \"ctDNA data var calling was done\" liubei\@cheerlandgroup.com\n";

sub mapping {
  my $sample = $_[0];
  return (
    "/usr/bin/bwa mem -M  -t 8 -R \"\@RG\\tID:$sample.sam\\tPL:ILLUMINA\\tSM:$sample.sam\\tDS:ref=$fa,pfx=$fa\" $fa $in/CleanData/$sample.R1.clean.fastq.gz $in/CleanData/$sample.R2.clean.fastq.gz| $samtools view -bS -t $fa.fai - > $out/bam/$sample.bam
$samtools fixmate -r -p $out/bam/$sample.bam $out/bam/$sample.fixmate.bam
$samtools sort -m 16G -@ 8 $out/bam/$sample.fixmate.bam -T $out/bam/$sample.sort -o $out/bam/$sample.sort.bam 
$samtools view -@ 8 -L $target -b -o $out/bam/$sample.target.bam $out/bam/$sample.sort.bam 
$samtools index $out/bam/$sample.sort.bam
$samtools index $out/bam/$sample.target.bam 
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../../bin/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -U ALLOW_SEQ_DICT_INCOMPATIBILITY  -R $fa -I $out/bam/$sample.target.bam -knownSites /data1/software/b37/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites /data1/software/b37/1000G_phase1.indels.b37.vcf -knownSites /data1/software/b37/dbsnp_138.b37.vcf -o $out/bam/$sample.target.recal.table -rf BadCigar -L /data1/workdir/liubei/ctDNA/20170503_11sample_test/cosmic.list.bed
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../../bin/GenomeAnalysisTK.jar -T PrintReads -nct 8 -R $fa -I $out/bam/$sample.target.bam -BQSR $out/bam/$sample.target.recal.table -o $out/bam/$sample.target.recal.bam -rf BadCigar
$samtools index $out/bam/$sample.target.recal.bam
"
  );
}

sub rmdup {
  my $sample = $_[0];
  my @chr = ( 1 .. 22, "X" );
  my $result;
  foreach my $chr ( 0 .. $#chr ) {
    $result
      .= "perl $bin/../index/index.count.pl -in $out/bam/$sample.target.recal.bam -out $out/bam/$sample.target.recal.index -chr $chr > $out/bam/$sample.target.recal.index.fiter2\n";
    $result
      .= "sort -k4n -k5n -k1 -k2 $out/bam/$sample.target.recal.index.chr$chr.log > $out/bam/$sample.target.recal.index.chr$chr.log.sort \n";
    $result
      .= "perl $bin/../index/index.filter.pl -in $out/bam/$sample.target.recal.index.chr$chr.log.sort -out $out/bam/$sample.target.recal.index.chr$chr.log.sort.rmdup\n";
  }
  return ("$result");
}

sub realign {
  my $sample = $_[0];
  return (
    "cat $out/bam/$sample.target.recal.index.chr*.log.sort.rmdup  > $out/bam/$sample.target.recal.index.reads
$samtools view -h $out/bam/$sample.target.recal.bam > $out/bam/$sample.target.recal.sam
perl $bin/../index/collect.pl -list $$out/bam/$sample.target.recal.index.reads -out $out/bam/$sample.target.recal.index.rmdup.sam -sam $out/bam/$sample.target.recal.sam
$samtools view -S $out/bam/$sample.target.recal.index.rmdup.sam -o $out/bam/$sample.target.recal.index.rmdup.bam -b
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../../bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -U ALLOW_SEQ_DICT_INCOMPATIBILITY -nt 8 -R $fa -I $out/bam/$sample.target.recal.index.rmdup.bam -known /data1/software/b37/Mills_and_1000G_gold_standard.indels.b37.vcf -known /data1/software/b37/1000G_phase1.indels.b37.vcf  -o $out/bam/$sample.target.intervals.list -rf BadCigar 
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../../bin/GenomeAnalysisTK.jar -T IndelRealigner -U ALLOW_SEQ_DICT_INCOMPATIBILITY -model USE_READS -R $fa -I $out/bam/$sample.target.recal.index.rmdup.bam -known /data1/software/b37/Mills_and_1000G_gold_standard.indels.b37.vcf -known /data1/software/b37/1000G_phase1.indels.b37.vcf -targetIntervals $out/bam/$sample.target.intervals.list -o $out/bam/$sample.target.realigned.bam -rf BadCigar
$samtools index $out/bam/$sample.target.realigned.bam
rm  $out/bam/$sample.target.recal.bam   
rm $out/bam/$sample.target.recal.sam $out/bam/$sample.target.recal.index.rmdup.sam  
"
  );
}

sub bam_stat {
  my $sample = $_[0];
  return (
    "$samtools flagstat $out/bam/$sample.sort.bam > $out/bam.stat/$sample.sort.bam.flagstat
$samtools flagstat $out/bam/$sample.target.bam > $out/bam.stat/$sample.target.bam.flagstat
$samtools depth -q 20 -Q 20 -l 60 -d 500000 $out/bam/$sample.target.bam >$out/DP/$sample.target.bam.depth
$samtools depth -q 20 -Q 20 -l 60 -d 500000 $out/bam/$sample.target.realigned.bam > >$out/DP/$sample.target.realigned.bam.depth
perl $bin/cov.pl $target $out/DP/$sample.target.bam.depth
ln -s $out/bam/$sample.target.realigned.removemultiple.bam $out/release/bam/.
ln -s $out/bam/$sample.target.realigned.removemultiple.bam.bai $out/release/bam/.
"
  );
}

sub cnv {
  my $sample = $_[0];
  return ("sh $bin/../cnv/cnv.sh $out");
}

sub snp_calling {
  my ( $sample, $tumor ) = @_;
  my @chr = ( 1 .. 22, "X" );
  my @return_vcf;
  open my $LIST, ">", "$out/vcf/$sample.list";
  foreach (@chr) {
    my $chr = $_;
    push @return_vcf,
      "java -Xmx10g -Djava.io.tmpdir=/tmp -jar $bin/../../bin/GenomeAnalysisTK.jar -nct 2 -T MuTect2 -R $fa -I:tumor  $out/bam/$sample.target.realigned.removemultiple.bam -o $out/vcf/$sample.chr$chr.MuTect2.vcf -hets 0.0001 -L $chr";
    print $LIST "$out/vcf/$sample.chr$chr.MuTect2.vcf\n";
  }
  my $return_vcf = join "\n", @return_vcf;
  return (
    "$return_vcf\nperl $bin/filter.bam.test.pl -in $out/bam/$sample.target.realigned.removemultiple.bam -out $out/bam/$sample.bam.filter"
  );
}

sub anno {
  my ( $sample, $tumor ) = @_;

  #  return("/usr/bin/java -Xmx40g -Djava.io.tmpdir=/tmp -jar /data1/software/ReSeq/GATK-3.7/GenomeAnalysisTK.jar -T MuTect2 -R $fa -I:tumor  $out/bam/$sample.target.realigned.removemultiple.bam -o $out/vcf/$sample.MuTect2.vcf -hets 0.0001
  return (
    "perl $bin/merge.vcf.pl -in $out/vcf/$sample.list -out $out/vcf/$sample.MuTect2.vcf
perl $bin/merge.cosmic.pl $out/vcf/$sample.MuTect2.vcf
perl $bin/../chemical_medicine/trans.vcf.chemical_medicine.pl $out/vcf/$sample.MuTect2.vcf $bin/../chemical_medicine/$tumor $out/result/$sample
perl $anno/convert2annovar.pl -format vcf4old $out/vcf/$sample.MuTect2.vcf  —includeinfo —withzyg —comment > $out/anno/$sample.annovar.txt
perl $anno/table_annovar.pl $out/anno/$sample.annovar.txt $anno/humandb -buildver hg19 -out $out/anno/$sample -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,clinvar_20160302,snp138,ljb26_all -operation g,r,r,f,f,f,f,f -nastring .
sed '1d' $out/anno/$sample.hg19_multianno.txt |cut -f 6-|paste -d \"\\t\" $out/anno/$sample.annovar.txt -> $out/anno/$sample.anno2.txt
perl $bin/merge.vcf.anno.pl $out/vcf/$sample.MuTect2.vcf.out $out/anno/$sample.anno2.txt
perl $bin/filter.somtic.pl $out/vcf/$sample.MuTect2.vcf.out $out/anno/$sample.anno2.txt.merge $out/result/$sample.somtic.xls
perl $bin/filter.bam.type.pl $out/bam/$sample.bam.filter $out/result/$sample.somtic.xls $out/result/$sample.somtic.new.xls $bin/../data_base/tran.id.new4
perl $bin/test.write.xlsx.pl $out/result/$sample.somtic.new.xls $out/result/$sample.somtic.xlsx
perl $bin/test.write.xlsx.pl $out/result/$sample.CHEMICAL.medicine.xls $out/result/$sample.CHEMICAL.medicine.xlsx
ln -s $out/result/$sample.somtic.xlsx $out/release/result/.
ln -s $out/result/$sample.CHEMICAL.medicine.xlsx $out/release/result/.
"
  );
}
