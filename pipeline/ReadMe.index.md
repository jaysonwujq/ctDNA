# ctDNA pipeline

for ctDNA analysis

## manual
```
perl pipeline/pip.index.pl 
-in inpath for cnv 
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
  PEmap -> DEunpair;
  Rmdup[label = "indexed data rmdup:\n\nfilter the reads have no right index\nget each index mapped to the same pos num\nget each type reads supports and the best reads name to respect the type\nget the type reads supports >=1,2,3,4 respectivly"];
  Rmdup[fillcolor=purple, style="rounded,filled", shape=box];
  Target [label = "get targeted bam"];
  DEunpair [label = "delete unmapped or \n unpaired reads"];
  STAT[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter1[fillcolor=yellow, style="rounded,filled", shape=box];
  Filter2[fillcolor=yellow, style="rounded,filled", shape=box];
  CHEMICAL_medicine_Type[fillcolor=yellow, style="rounded,filled", shape=diamond, shape=box];
  RmXA [label = "rm Multiple alignment"];
  SNVcalling [label = "calling snv by mutect2"];
  Add [label ="add cosmic and rs data"];
  STAT [label = "stat mismatch had high quality:\nmapping quality >=30\ninsersion and deletion < 2\nall type mismatch <= 3\nborder > 5bp\nbase quality >=20 or 80% nearby 10bp base quality >= 20"];
  Filter1 [label = "filter depth low:\ndepth all one base >= 30 \ndepth of alt >=5\nalt ratio >= 0.002"];
  Filter2 [label = "get mismatch had high quality:\n1.save the snv supported by >= 5 differ type reads\n and >= 1 positive read >= 1 reverse read\n2.indel have no high quality reads support tagged as non-filter\n3.snp have no high quality reads support tagged as filtered"];
  subgraph cluster0 {
    label = "mapping bwa";
    DEunpair -> Sorted;
    Sorted -> Target;
    Target -> Recal;
  }
  Recal -> RmXA;
  subgraph cluster1 {
    label = "bam filter";
    RmXA -> Rmdup;
  }
  Rmdup -> SNVcalling;
  subgraph cluster2{
    label = "SNV calling";
    SNVcalling -> Add
    Add -> ANNO;
  }
  subgraph cluster3{
    label = "somatic SNV filter";
    STAT -> Filter1;
    Filter1 -> Filter2;
    Filter2 -> "sample.somtic.new.xlsx";
  }
  splines=ortho;
  CHEMICAL_medicine_Type [label = "CHEMICAL medicine Type:\nhom_alt alt_ratio >= 0.8\nhet alt_ratio >= 0.3\nother hom_ref\nthe ssr is not detect,\n one of it is typed by indel"];
  ANNO -> STAT;
  Chemical -> CHEMICAL_medicine_Type;
  SNVcalling -> CHEMICAL_medicine_Type;
  CHEMICAL_medicine_Type -> "sample.CHEMICAL.medicine.xlsx";
  "sample.somtic.new.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  "sample.CHEMICAL.medicine.xlsx"[fillcolor=green, style="rounded,filled", shape=box];
  { rank=same;  "sample.somtic.new.xlsx","sample.CHEMICAL.medicine.xlsx";}
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
4. get targeted bam 
5. recal

#### filter 
for snp, indel, cnv calling
1. rm Multiple alignment
   * first, get the reads had tag XA
   * second, get both side mismatch bases number
   * third, classification types as mismatch 
      - 0 mismatch -> 1;
      - 1,2 mismatch -> 2;
      - >=3 mismacth -> 3;
      - soft-slip map +1;
   * last, delete best and second masth both the side have the same tpe
2. no md or cigar
3. mapping quality < 30
4. filter much mismatch
   * insersion and deletion >= 3
   * all type mismatch > 5
   * #soft clip and border 5bp uncount

#### rmdup
1. filter the reads have no right index
2. get each index mapped to the same pos num
3. get each type reads supports and the best reads name to respect the type
4. get the type reads supports >=1,2,3,4 respectivly

#### realign
1. get the rmdup bam
2. realign

#### bam stat
1. mapping reads 
2. target ratio
3. depth  target depth
4. coverage

### SNV

#### SNV  calling
1. calling snv by mutect2 

#### ANNO
1. add cosmic and rs data
2. annovar

#### somatic snv filter
1. stat mismatch had high quality 
   * mapping quality >=30
   * insersion and deletion < 2
   * all type mismatch <= 3
   * border > 5bp
   * base quality >=20 or 80% nearby 10bp base quality >= 20
2. filter depth low
   * depth all one base >= 30
   * depth of alt >=5
   * alt ratio >= 0.002
3. filter by high quality reads 
   * get mismatch had high quality
   * save the snv supported by >= 5 differ type reads and >= 1 positive read >= 1 reverse read
   * indel have no high quality reads support tagged as non-filter
   * snp have no high quality reads support tagged as filtered

#### CHEMICAL medicine
1. hom_alt alt_ratio >= 0.8
2. het alt_ratio >= 0.3
3. other hom_ref
4. the ssr is not detect, one of it is typed by indel

## dependents
   * bwa
   * samtools
   * GATK
   * Rscript
   * java

## need developments

