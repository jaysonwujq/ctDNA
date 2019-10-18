## for new chemical_medicine
  you should do 


```
  sh do.sh
```


# 1.from excel to text file for each type:lung ,all,jzc

example:
```
perl trans_from_excel_to_text.pl -in 化疗药物数据库更新-化药单药版-20171128.xlsx
```

# 2.do :
- perl tiqu.rs.pl  for above three file
- get the postion message for each vartation.


  example: 
  ```
  perl tiqu.rs.pl all 
  ```


# 3.do :
- perl extract_target.rs.pl  for above three file 
- get the targeted vars and trans the MSI to indels


example:
    ```
    perl extract_target.rs.pl all.RS
    ```
