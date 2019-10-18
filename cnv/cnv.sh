in=$1
cd $in/amplication/
bin=$2
cat $bin/control.list $in/list > $in/amplication/list

cat $bin/control.list |while read a;do
  ln -sf $bin/DP-control/$a.target.bam.depth $in/amplication/.
done

cat $in/list|while read a;do
  ln -sf $in/DP/$a.target.bam.depth $in/amplication/.
done

perl $bin/merge.pl $in/amplication/list $in/amplication/

cut -f 5 $bin/cosmic.list.bed.old2.OUT.sort2 |uniq|while read a;do
  Rscript $bin/ctDNA.R $a $in/amplication/
done
rm -f $in/amplication/list.check
touch $in/amplication/list.check
cat $in/amplication/list|while read a;do
  cut -f 5 $bin/cosmic.list.bed.old2.OUT.sort2 |uniq|while read b;do
    echo "$a.$b" >>$in/amplication/list.check  
  done
done
rm -f  $in/cnv.out
perl $bin/filter.cnv.pl $in/amplication/list.check $in/amplication $in/cnv.out
