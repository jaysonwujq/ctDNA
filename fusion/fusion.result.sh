bin=$1
in=$2
ln -sf $bin/../fusion/batch6_dd-42_S4.predSV.txt.filter.ANNO.filtered $in/.
grep . $in/*ANNO.filtered|grep -v type |sed -e 's/.predSV.txt.filter.ANNO.filtered:/\t/g' |awk -F '/' '{print $NF}' > $in/fusion.result
cat $bin/../fusion/head.new $in/fusion.result > $in/fusion.all.result
perl $bin/test.write.xlsx.pl $in/fusion.all.result $in/../release/result/fusion.result.xlsx
