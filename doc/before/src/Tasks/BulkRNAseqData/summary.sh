#!/bin/bash

result_path=${1}
file_gtf=${2}
file_current=`realpath ${BASH_SOURCE}`
dir_file_current=`dirname ${file_current}`


get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

# 3.Expression quantification ----------------------------------------------------------------------
echo "3.Expression quantification..."

echo "121 `get_time` RNA_seq_preprocessing ..."
/usr/bin/python3.8 ${dir_file_current}/RNA_seq_preprocessing.py ${result_path}

echo "124 `get_time` RemoveDuplicateGenes ..."
/usr/bin/python3.8 ${dir_file_current}/RemoveDuplicateGenes.py ${result_path}/merge_fpkm.txt ${result_path}/GeneSymbol_fpkm.txt
mv ${result_path}/merge_fpkm.txt ${result_path}/Transcript_fpkm.txt

echo "134 `get_time` prepDE ..."
/usr/bin/python2 ${dir_file_current}/prepDE.py -i ${result_path}/gtf_file_path.txt -g ${result_path}/gene_count.txt -t ${result_path}/transcript_count.txt > ${result_path}/log/CountConvertion.log 2>&1
sed -i "s/,/\t/g" ${result_path}/gene_count.txt
sed -i "s/,/\t/g" ${result_path}/transcript_count.txt
mv ${result_path}/gene_count.txt ${result_path}/GeneSymbol_counts.txt
mv ${result_path}/transcript_count.txt ${result_path}/Transcript_counts.txt


rm -rf ${result_path}/GeneSymbol_counts.txt
rm -rf ${result_path}/GeneSymbol_fpkm.txt
/usr/bin/python3.8 ${dir_file_current}/gfold_fpkm.py ${result_path}/bam ${result_path} ${file_gtf}