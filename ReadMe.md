# ctDNA pipeline![]()

for ctDNA analysis

## manual
```
perl pipeline/pip.clean.pl 
-in inpath data 
-out outpath 
-target target bed (defult $bin/../data_base/first.panel.bed)
-fa fa (defult /data1/software/b37/human_g1k_v37.fasta)
-db dataBase (defult /data1/software/b37)
-anno anno database (defult /data1/software/annovar/)
```
## pipeline description

flowchart
```plantuml
digraph G {
  compound=true;
  CleanData -> Chemical;
  Chemical [label = "get CHEMICAL medicine type"];
  CleanData -> PEmap;
  PEmap [label = "paired reads mapping"];
  raw->mapping->bam_filter->var_calling->var_filter->result[arrowhead=none];
  PEmap -> DEunpair;
  Rmdup1 [label = "rmdup for fusion"];
  Realign [label = "recal && realign"];
  Target [label = "get targeted bam"];
  DEunpair [label = "delete unmapped or \n unpaired reads"];
  FilterOther[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter[fillcolor=yellow, style="rounded,filled", shape=box];
  STAT[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter1[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter2[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter_fusion1[fillcolor=yellow, style="rounded,filled", shape=box];
  ANNO_fusion[shape=box];
  CHEMICAL_medicine_Type[fillcolor=yellow, style="rounded,filled", shape=diamond, shape=box];
  RmXA [label = "rm Multiple alignment"];
  FilterOther [label ="filter no md or cigar\nfilter much mismatch:\ninsersion and deletion >= 3\nall type mismatch > 5\nsoft clip and border 5bp uncount"];
  MergeDP [label = "merge business and control depth"];
  Filter [label = "filter cnv by gene:\nmore than 70% base of one gene > 0.7\nbase of one gene detected base > 500"];
  Chrcalling [label = "calling each chr snv by mutect2"];
  Merge [label = "merge vcf"];
  Add [label ="add cosmic and rs data"];
  STAT [label = "stat mismatch had high quality:\nmapping quality >=30\ninsersion and deletion < 2\nall type mismatch <= 3\nborder > 5bp\nbase quality >=20 or 80% nearby 10bp base quality >= 20"];
  Filter1 [label = "filter depth low:\ndepth all one base >= 30 \ndepth of alt >=5\nalt ratio >= 0.002"];
  Filter2 [label = "get mismatch had high quality:\n1.save the snv supported by >= 5 differ type reads\n and >= 1 positive read >= 1 reverse read\n2.indel have no high quality reads support tagged as non-filter\n3.snp have no high quality reads support tagged as filtered"];
  FUSION [label = "CREST detect SV raw result"];
  Filter_fusion1 [label = "filter:\none or more percent_identity < 0.9\none or more depth < 4\nboth side out of target +-50bp"];
  POS [label = "get the postion of both split side"];
  ANNO_fusion [label = "annovar for the both side split pos\n to get gene information"];
  merge_result [label = "merge the annovar result to sv result"];
  subgraph cluster0 {
    label = "mapping bwa";
    DEunpair -> Sorted;
    Sorted -> Target;
    Target -> Realign;
    Sorted -> Rmdup1;
    Realign -> Rmdup;
    Rmdup1 -> "sample.rmdup.bam"[splines="line"];
    Realign -> "sample.target.realigned.rmxa.bam";
    {rank=same;  Rmdup, Realign};
  }
  subgraph cluster1 {
    label = "bam filter";
    RmXA -> FilterOther;
    FilterOther -> "sample.target.realigned.removemultiple.bam";
  }
  subgraph cluster3{
    label = "CNV calling";
    "sample.target.bam.depth" -> MergeDP;
    MergeDP -> CNVcalling;
    CNVcalling -> Filter;
  }
  subgraph cluster4{
    label = "SNV calling";
    Chrcalling -> Merge;
    Merge -> Add;
    Add -> ANNO;
  }
  subgraph cluster5{
    label = "somatic SNV filter";
    STAT -> Filter1;
    Filter1 -> Filter2;
    Filter2 -> "sample.somtic.new.xlsx";
  }
  subgraph cluster6{
    label = "FUSION detect by CREST";
    FUSION -> Filter_fusion1;
    Filter_fusion1 -> POS;
    POS -> ANNO_fusion;
    ANNO_fusion -> merge_result;
  }
  splines=ortho;
  CHEMICAL_medicine_Type [label = "CHEMICAL medicine Type:\nhom_alt alt_ratio >= 0.8\nhet alt_ratio >= 0.3\nother hom_ref\nthe ssr is not detect,\n one of it is typed by indel"];
  "sample.rmdup.bam" -> FUSION;
  "sample.target.realigned.rmxa.bam" -> RmXA;
  "sample.target.realigned.removemultiple.bam" -> Chrcalling;
  Target -> "sample.target.bam.depth";
  ANNO -> STAT;
  Chemical -> CHEMICAL_medicine_Type;
  Merge -> CHEMICAL_medicine_Type;
  CHEMICAL_medicine_Type -> "sample.CHEMICAL.medicine.xlsx";
  merge_result -> "fusion.result.xlsx";
  Filter -> "cnv.xlsx";
  "cnv.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  "sample.somtic.new.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  "sample.CHEMICAL.medicine.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  "fusion.result.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  "sample.target.realigned.rmxa.bam"[fillcolor=pink, style="filled",shape=box];
  "sample.target.realigned.removemultiple.bam"[fillcolor=pink, style="filled",shape=box];
  "sample.target.bam.depth"[fillcolor=pink, style="filled",shape=box];
  "sample.rmdup.bam"[fillcolor=pink, style="filled",shape=box];
  { rank=same; raw; CleanData;}
  { rank=same; bam_filter; "sample.target.realigned.rmxa.bam", "sample.target.bam.depth";}
  { rank=same; var_calling; "sample.rmdup.bam","sample.target.realigned.removemultiple.bam";}
  { rank=same; var_filter;  FUSION;}
  { rank=same; result; "sample.somtic.new.xlsx","sample.CHEMICAL.medicine.xlsx","fusion.result.xlsx","cnv.xlsx";}
}
```

      
### get CHEMICAL medicine type
* 1 for 肿瘤个体化用药
* 2 for 肺癌
* 3 for 结直肠癌

### mapping
#### bwa
1. paired reads mapping
2. delete unmapped reads and unpaired reads
3. sort bam
4. get targeted bam ; for snp, indel, cnv calling  
5. recal
6. realign
7. rmdup 
   * after 6; __just stat duplication__
   * after 3; for fusion detect 

#### filter 
for snp, indel, cnv calling
1. rm Multiple alignment
2. no md or cigar
3. mapping quality < 30
4. filter much mismatch(__too much mismatch in one read maybe caused by wrong mapping__)
   * insersion and deletion >= 3
   * all type mismatch > 5
   * #soft clip and border 5bp uncount

#### bam stat
1. mapping reads 
2. target ratio (__evalulate the panel Capture efficiency__)
3. depth  target depth 
4. coverage (__differ depth coverage can reflect can it be used to calling snv__)

### CNV
#### CNV calling 
1. get depth of target region
2. merge business and control depth 
3. cnv calling for each base

#### CNV result 
1. filter 
   * mosre than 70% base of one gene > 0.7 (__detect gene as unit__)
   * bases of one gene > 500 (__too little have too much false positives__)

### SNV

#### SNV  calling
1. calling each chr snv by mutect2 (__because it can both save time and avoid GATK bug when it do much time__)
2. merge vcf 

#### ANNO
1. add cosmic and rs data
2. annovar

#### somatic snv filter
1. stat mismatch had high quality (__prepare filter false positives__)
   * mapping quality >=30
   * insersion and deletion < 2
   * all type mismatch <= 3
   * border > 5bp
   * base quality >=20 or 80% nearby 10bp base quality >= 20

2. filter depth low (__low depth or ratio generaly means false positives__)
   * depth all one base >= 30
   * depth of alt >=5
   * alt ratio >= 0.002
3. filter by high quality reads (__because not do rmdup, to avoid snv origin is PCR__)
   * get mismatch had high quality
   * save the snv supported by >= 5 differ type reads and >= 1 positive read >= 1 reverse read
   * indel have no high quality reads support tagged as non-filter ( __because mutect2 do cluster and Stitching, the postion may changed__)
   * snp have no high quality reads support tagged as filtered

#### CHEMICAL medicine 
1. hom_alt alt_ratio >= 0.8 (__homozygous mutation__) 
2. het 0.8 > alt_ratio >= 0.3 (__heterozygote mutation__)
3. other hom_ref (__ref type homozygous__)
4. the ssr is not detect, one of it is typed by indel

### FUSION
1. SV detect by CREST
2. filter 
   * one or more percent_identity < 0.9 (__maybe wrong mapped postion__)
   * one or more depth < 4 (__maybe wrong or ratio is very low__)
   * both side out of target +-50bp (__one side must in the target, start and end expand 50bp__)
3. anno (__the anno software can not anno varation between different chr__)
   * get the postion of both side
   * annovar for the both side 
   * merge the annovar result to sv result
4. anno fiter
   * filtered-intronic-same(__just SV, not gene fusion__)
   * filtered-low_coverage both side < 300 (__one side must in the target, start and end expand 50bp and depth >= 300__ )
   * filtered-intergenic(__just SV, not gene fusion__)
   * filtered-strand(__just SV, not gene fusion__)

## dependents
   * bwa
   * samtools
   * GATK
   * CREST
   * Rscript
   * java
   * blat server (/gfServer start 192.168.1.205 8000 /home/liubei/bin/x86_64/CREST/human_g1k_v37.fasta.2bit)

## need developments
1. some data have some problems for experimental stage may can not done or need very long time for crest
2. cnv can splited to two type plasma and others 
3. SSR detect ?

## dependents perl packages
   * Bio::DB::Sam;
   * Data::Dumper;
   * Excel::Writer::XLSX;
   * File::Basename
   * File::Spec;
   * FindBin
   * Getopt::Long;
   * IO::All
   * List::Util
   * Spreadsheet::XLSX;
   * Statistics::Basic
   * Text::Iconv;

## dependents R library
   * cghFLasso
