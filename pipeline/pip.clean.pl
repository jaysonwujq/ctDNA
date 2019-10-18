#!/bin/usr/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ( $in, $out, $list, $target, $var, $type, $fa, $anno, $data_base );
GetOptions(
  "in:s"     => \$in,
  "out:s"    => \$out,
  "target:s" => \$target,
  "fa:s"     => \$fa,
  "anno:s"   => \$anno,
  "db:s"     => \$data_base,
);
my $bin      = $Bin;
my $samtools = "$bin/../bin/samtools";

`mkdir -p $out $out/shell $out/clean $out/bam $out/vcf $out/stat $out/bam.stat $out/anno $out/amplication $out/fusion $out/result $out/DP`;
`mkdir -p $out/release $out/release/bam $out/release/result`;
`perl $bin/from_sample_to_type.pl -in $in -out $out`;
`ln -sf $in/stats/list $out/list`;
$type = "$out/type";
if ( !$fa )   { $fa   = "/data1/software/b37/human_g1k_v37.fasta"; }
if ( !$data_base )   { $data_base   = "/data1/software/b37"; }
if ( !$anno ) { $anno = "/data1/software/annovar/"; }
if ( !$target){ $target = "$bin/../data_base/first.panel.bed";}
open my $SNP_SH, ">", "$out/shell/raw_to_vcf.sh";
open my $DO_SH,  ">", "$out/shell/do.qsub.sh";

open my $LIST, "<", "$in/stats/list";
open my $TYPE, "<", "$type";
my $i = 0;
my $tumor_type;
my %case;
my %control;
my %sample;

my $count = 0;
while (<$LIST>) {
  $count++;
  chomp;
  open my $SH1,       ">", "$out/shell/$_.sh";
  open my $SHBWA,     ">", "$out/shell/$_.bwa.sh";
  open my $DOBWA,     ">", "$out/shell/do.$_.bwa.sh";
  open my $SHBWASTAT, ">", "$out/shell/$_.bwa.stat.sh";
  open my $SNP,       ">", "$out/shell/$_.snpcalling.sh";
  open my $ANNO,      ">", "$out/shell/$_.anno.sh";
  open my $FUSION,    ">", "$out/shell/$_.fusion.sh";
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
  my $bwa = &mapping($_);
  print $SHBWA "$bwa\n";
  my $bwa_stat    = &bam_stat($_);
  my $snp_calling = &snp_calling( $_, $tumor_type );
  my $anno_result = &anno( $_, $tumor_type );
  my $fusion      = &fusion ($_);
  print $FUSION "$fusion";
  print $SHBWASTAT "$bwa_stat\nsh $bin/../cnv/cnv.sh $out $bin/../cnv\nperl $bin/test.write.xlsx.pl $out/cnv.out $out/release/result/cnv.xlsx\n";
  print $SNP "$snp_calling\nsh $out/shell/$_.bwa.stat.sh\n";
  print $ANNO "$anno_result\n";
  print $ANNO "sh $out/shell/do.qsub.sh\n";
  print $DOBWA "sh $out/shell/$_.bwa.sh\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=8,mem=90G,walltime=1000:00:00 --maxjob 50 $out/shell/do.$_.bwa.sh --convert n --jobprefix ct.$i.bwa --queue batch \n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 1  --resource nodes=1:ppn=2,mem=10G,walltime=1000:00:00 --maxjob 50 $out/shell/$_.snpcalling.sh --convert n --jobprefix ct.$i.snp --queue all\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=1,mem=10G,walltime=1000:00:00 --maxjob 50 $out/shell/$_.anno.sh --convert n  --jobprefix ct.$i.anno --queue batch\n";
  print $SH1
    "perl $bin/qsub-pbs.pl --lines 50 --resource nodes=1:ppn=1,mem=10G,walltime=10:00:00 --maxjob 50 $out/shell/$_.fusion.sh --convert n  --jobprefix ct.$i.fusion --queue batch\n";
  close $SH1;
  print $SNP_SH "sh $out/shell/$_.sh&\n";
}

print $DO_SH
  "perl $bin/target.dp.pl $in/stats/list $target $out/DP/all.dp.stat $out/DP\n";
print $DO_SH "perl $bin/stat.target.pl $in/stats/list $out/bam.stat\n";
print $DO_SH
  "perl $bin/merge.all.stat.pl $in/stats/clean.stat $out/bam.stat/target.ratio $out/DP/all.dp.stat $out/all.stat\n";
print $DO_SH
  "perl $bin/test.write.xlsx.pl $out/all.stat $out/release/result/all.stat.xlsx\n";
`sh $out/shell/raw_to_vcf.sh`;
sub mapping {
  my $sample = $_[0];
  return (
    "/usr/bin/bwa mem -M  -t 8 -R \"\@RG\\tID:$sample.sam\\tPL:ILLUMINA\\tSM:$sample.sam\\tDS:ref=$fa,pfx=$fa\" $fa $in/CleanData/$sample.R1.clean.fastq.gz $in/CleanData/$sample.R2.clean.fastq.gz| $samtools view -bS -t $fa.fai - > $out/bam/$sample.bam
$samtools fixmate -r -p $out/bam/$sample.bam $out/bam/$sample.fixmate.bam
$samtools sort -m 16G -@ 8 $out/bam/$sample.fixmate.bam -T $out/bam/$sample.sort -o $out/bam/$sample.sort.bam 
$samtools view -@ 8 -L $target -b -o $out/bam/$sample.target.bam $out/bam/$sample.sort.bam 
$samtools index $out/bam/$sample.sort.bam
$samtools index $out/bam/$sample.target.bam 
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../bin/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -U ALLOW_SEQ_DICT_INCOMPATIBILITY  -R $fa -I $out/bam/$sample.target.bam -knownSites $data_base/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites $data_base/1000G_phase1.indels.b37.vcf -knownSites $data_base/dbsnp_138.b37.vcf -o $out/bam/$sample.target.recal.table -rf BadCigar -L $target.change.mode.bed
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../bin/GenomeAnalysisTK.jar -T PrintReads -nct 8 -R $fa -I $out/bam/$sample.target.bam -BQSR $out/bam/$sample.target.recal.table -o $out/bam/$sample.target.recal.bam -rf BadCigar
$samtools index $out/bam/$sample.target.recal.bam
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -U ALLOW_SEQ_DICT_INCOMPATIBILITY -nt 8 -R $fa -I $out/bam/$sample.target.recal.bam -known $data_base/Mills_and_1000G_gold_standard.indels.b37.vcf -known $data_base/1000G_phase1.indels.b37.vcf  -o $out/bam/$sample.target.intervals.list -rf BadCigar 
/usr/bin/java -Xmx20g -Djava.io.tmpdir=/tmp -jar $bin/../bin/GenomeAnalysisTK.jar -T IndelRealigner -U ALLOW_SEQ_DICT_INCOMPATIBILITY -model USE_READS -R $fa -I $out/bam/$sample.target.recal.bam -known $data_base/Mills_and_1000G_gold_standard.indels.b37.vcf -known $data_base/1000G_phase1.indels.b37.vcf -targetIntervals $out/bam/$sample.target.intervals.list -o $out/bam/$sample.target.realigned.bam -rf BadCigar
$samtools index $out/bam/$sample.target.realigned.bam
java -Xmx20g -jar $bin/../bin/MarkDuplicates.jar  INPUT=$out/bam/$sample.target.realigned.bam OUTPUT=$out/bam/$sample.target.realigned.rmdup.bam METRICS_FILE=$out/bam/$sample.target.realigned.rmdup.metrix  CREATE_INDEX=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
java -Xmx20g -jar $bin/../bin/MarkDuplicates.jar  INPUT=$out/bam/$sample.sort.bam OUTPUT=$out/bam/$sample.rmdup.bam METRICS_FILE=$out/bam/$sample.rmdup.metrix  CREATE_INDEX=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
samtools index $out/bam/$sample.target.realigned.rmdup.bam
samtools index $out/bam/$sample.rmdup.bam
$samtools view -@ 8 -h $out/bam/$sample.target.realigned.bam > $out/bam/$sample.target.realigned.sam
perl $bin/filter.bam.xa.pl -i $out/bam/$sample.target.realigned.bam -bam $out/bam/$sample.target.realigned.sam -o $out/bam/$sample.target.realigned.rmxa.sam
$samtools view -@ 8 -S $out/bam/$sample.target.realigned.rmxa.sam -o $out/bam/$sample.target.realigned.rmxa.bam -b
$samtools index $out/bam/$sample.target.realigned.rmxa.bam
perl $bin/filter.bam.pl -i $out/bam/$sample.target.realigned.rmxa.bam -bam $out/bam/$sample.target.realigned.rmxa.sam -o $out/bam/$sample.target.realigned.filter.sam
$samtools view -@ 8 -S $out/bam/$sample.target.realigned.filter.sam -o $out/bam/$sample.target.realigned.removemultiple.bam -b
$samtools index $out/bam/$sample.target.realigned.removemultiple.bam
rm  $out/bam/$sample.target.recal.bam $out/bam/$sample.target.realigned.bam  
rm $out/bam/$sample.target.realigned.filter.sam  $out/bam/$sample.target.realigned.sam
"
  );
}

sub bam_stat {
  my $sample = $_[0];
  return (
    "$samtools flagstat $out/bam/$sample.sort.bam > $out/bam.stat/$sample.sort.bam.flagstat
$samtools flagstat $out/bam/$sample.target.bam > $out/bam.stat/$sample.target.bam.flagstat
$samtools depth -q 20 -Q 20 -l 60 -d 500000 $out/bam/$sample.target.bam >$out/DP/$sample.target.bam.depth
perl $bin/cov.pl $target $out/DP/$sample.target.bam.depth
ln -s $out/bam/$sample.target.realigned.removemultiple.bam $out/release/bam/.
ln -s $out/bam/$sample.target.realigned.removemultiple.bam.bai $out/release/bam/.
#perl $bin/stat.single.target.pl $sample $out/bam.stat;
"
  );
}

sub cnv {
  my $sample = $_[0];
  return ("sh $bin/../cnv/cnv.sh $out $bin/../cnv $bin/../cnv");
}

sub snp_calling {
  my ( $sample, $tumor ) = @_;
  my @chr = ( 1 .. 22, "X" );
  my @return_vcf;
  open my $LIST, ">", "$out/vcf/$sample.list";
  foreach (@chr) {
    my $chr = $_;
    push @return_vcf,
      "java -Xmx10g -Djava.io.tmpdir=/tmp -jar $bin/../bin/GenomeAnalysisTK.jar -nct 2 -T MuTect2 -R $fa -I:tumor  $out/bam/$sample.target.realigned.removemultiple.bam -o $out/vcf/$sample.chr$chr.MuTect2.vcf -hets 0.0001 -L $chr";
    print $LIST "$out/vcf/$sample.chr$chr.MuTect2.vcf\n";
  }
  my $return_vcf = join "\n", @return_vcf;
  return (
    "$return_vcf\nperl $bin/filter.bam.test.pl -in $out/bam/$sample.target.realigned.removemultiple.bam -out $out/bam/$sample.bam.filter"
  );
}

sub anno {
  my ( $sample, $tumor ) = @_;
  return (
    "perl $bin/merge.vcf.pl -in $out/vcf/$sample.list -out $out/vcf/$sample.MuTect2.vcf
perl $bin/merge.cosmic.pl $out/vcf/$sample.MuTect2.vcf
perl $bin/../chemical_medicine/trans.vcf.chemical_medicine.pl $out/vcf/$sample.MuTect2.vcf $bin/../chemical_medicine/$tumor $out/result/$sample
perl $anno/convert2annovar.pl -format vcf4old $out/vcf/$sample.MuTect2.vcf  —includeinfo —withzyg —comment > $out/anno/$sample.annovar.txt
perl $anno/table_annovar.pl $out/anno/$sample.annovar.txt $anno/humandb -buildver hg19 -out $out/anno/$sample -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,clinvar_20160302,snp138,ljb26_all -operation g,r,r,f,f,f,f,f -nastring .
sed '1d' $out/anno/$sample.hg19_multianno.txt |cut -f 6-|paste $out/anno/$sample.annovar.txt -> $out/anno/$sample.anno2.txt
perl $bin/merge.vcf.anno.pl $out/vcf/$sample.MuTect2.vcf.out $out/anno/$sample.anno2.txt
perl $bin/filter.somtic.pl $out/vcf/$sample.MuTect2.vcf.out $out/anno/$sample.anno2.txt.merge $out/result/$sample.somtic.xls
perl $bin/filter.bam.type.pl -in $out -sample $sample -target $target
perl $bin/test.write.xlsx.pl $out/result/$sample.somtic.new.xls $out/release/result/$sample.somtic.xlsx
perl $bin/test.write.xlsx.pl $out/result/$sample.CHEMICAL.medicine.xls $out/release/result/$sample.CHEMICAL.medicine.xlsx
"
  );
}

sub fusion {
  my ( $sample ) = $_[0];
  return (
"perl /data1/workdir/liubei/ctDNA/crest/CREST/extractSClip.pl -i $out/bam/$sample.rmdup.bam --ref_genome $fa -o $out/fusion
perl /home/liubei/bin/x86_64/CREST/CREST.pl -f  $out/fusion/$sample.rmdup.bam.cover -d  $out/bam/$sample.rmdup.bam --ref_genome /data1/workdir/liubei/ctDNA/crest/CREST/human_g1k_v37.fasta -t /data1/workdir/liubei/ctDNA/crest/CREST/human_g1k_v37.fasta.2bit -o $out/fusion -p $sample 
perl $bin/../fusion/change.mode.for.anno.pl $out/fusion/$sample.predSV.txt $target
perl $anno/annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile $out/fusion/$sample -exonsort $out/fusion/$sample.predSV.txt.filter.anno.txt /data1/software/annovar//humandb
perl $bin/../fusion/anno.merge.pl $out/fusion/$sample.predSV.txt.filter $out/fusion/$sample.variant_function
perl $bin/../fusion/merge.exon.intro.pl $out/fusion/$sample.predSV.txt.filter.ANNO ~/workdir/program/cds/genes.gtf.long.trns.exon2
sh $bin/../fusion/fusion.result.sh $bin $out/fusion/
"
  );
}
