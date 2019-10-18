excel=化疗药物数据库更新-化药单药版-20171128.xlsx
in=/home/liubei/workdir/pipline/program/ctDNA/chemical_medicine

perl $in/trans_from_excel_to_text.pl -in $in/$excel
perl $in/tiqu.rs.pl all 
perl $in/tiqu.rs.pl lung
perl $in/tiqu.rs.pl jzc
perl $in/extract_target.rs.pl all.RS
perl $in/extract_target.rs.pl lung.RS
perl $in/extract_target.rs.pl jzc.RS
