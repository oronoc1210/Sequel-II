cat */*_stats.txt | head -n 1 > Summarized_Results.xls
while read L
do
cat  */$L\_stats.txt |grep -v "subread" >> Summarized_Results.xls
done < libs.txt
