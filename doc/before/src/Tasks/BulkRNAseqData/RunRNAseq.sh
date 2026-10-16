#!/bin/bash

# Description : RNA-seq data analysis pipeline.
# Author      : Chengchen Zhao
# Version     : 1.1
# Time        : 2021-07-21 09:00:40

####################################################################################################
#################################     PARAMETER/FUNCTION DEFINE     ################################
####################################################################################################
data_path=${1}
result_path=${2}
out_name=${3}

layout=${4} # single/pair
fastq_file=${5} # "YourFastqFile.fastq.gz" for single-end
                # "YourFastqFile_R1.fastq.gz,YourFastqFile_R2.fastq.gz" for pair-end

mapping_index=${6} # /data/zhaochengchen/Data/mm9/mm9
refseq_gtf=${7} # /data/zhaochengchen/Data/mm9/mm9_refSeq.gtf
refseq_bed=${8} # /data/zhaochengchen/Data/mm9/mm9.refseq.bed
chromInfo=${9} # /data/zhaochengchen/Data/mm9/mm9.len

file_current=`realpath ${BASH_SOURCE}`
dir_file_current=`dirname ${file_current}`

task_num=$cdesk_task_num
if [ ! $cdesk_task_num ]; then
    task_num=8
fi

echo "--------------------------------------------------------------------------------" 
if [ ! -n "${9}" ]
then 
    echo "Parameters are missing, please check your input." 
    exit 2
else 
    echo "RNA-seq data analysis pipeline is now running...  ${task_num} ${data_path}  ${result_path}  ${out_name}  ${layout}  ${fastq_file}  ${mapping_index}  ${refseq_gtf}  ${refseq_bed}  ${chromInfo}"
fi 


log_path=${result_path}/log
bam_path=${result_path}/bam

####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################

if [ ! -d "${result_path}" ]; then
  mkdir ${result_path}
fi
if [ ! -d "${log_path}" ]; then
  mkdir ${log_path}
fi
if [ ! -d "${bam_path}" ]; then
  mkdir ${bam_path}
fi
if [ ! -d "${result_path}/${out_name}" ]; then
  mkdir ${result_path}/${out_name}
fi
if [ ! -d "${log_path}/${out_name}" ]; then
  mkdir ${log_path}/${out_name}
fi


get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}


# 1.mapping ----------------------------------------------------------------------------------------
echo "1.mapping..."
if [ "${layout}" == "single" ]
then
    echo "67 `get_time` trim_galore ..."
    trim_galore --gzip -j ${task_num} ${data_path}/${fastq_file} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    
    echo "70 `get_time` hisat2 ..."
    hisat2 -p ${task_num} --dta -x ${mapping_index} -U ${result_path}/${out_name}_trimmed.fq.gz -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${result_path}/${out_name}_trimmed.fq.gz
    
    echo "74 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam
    
    echo "77 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${result_path}/${out_name}/${out_name}.Sorted.bam
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam
    mv ${result_path}/${out_name}/${out_name}.Sorted.bam ${bam_path}/${out_name}.bam
    
    echo "82 `get_time` stringtie ..."
    stringtie -e -B -p ${task_num} -G ${refseq_gtf} -o ${result_path}/${out_name}/${out_name}.gtf ${bam_path}/${out_name}.bam > ${log_path}/${out_name}/${out_name}.Stringtie.log 2>&1
elif [ "${layout}" == "pair" ]
then
    echo "86 `get_time` trim_galore ..."
    oldIFS=$IFS
    IFS=','
    read -ra fastq_file_arr <<< "${fastq_file}"
    trim_galore --gzip -j ${task_num} --paired ${data_path}/${fastq_file_arr[0]} ${data_path}/${fastq_file_arr[1]} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    IFS=$oldIFS
    
    echo "91 `get_time` hisat2 ..."
    hisat2 -p ${task_num} --dta -x ${mapping_index} -1 ${result_path}/${out_name}_val_1.fq.gz -2 ${result_path}/${out_name}_val_2.fq.gz -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${result_path}/${out_name}_val_1.fq.gz ${result_path}/${out_name}_val_2.fq.gz
    
    echo "95 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam
    
    echo "98 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${result_path}/${out_name}/${out_name}.Sorted.bam
    
    echo "101 `get_time` stringtie ..."
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam
    mv ${result_path}/${out_name}/${out_name}.Sorted.bam ${bam_path}/${out_name}.bam
    stringtie -e -B -p ${task_num} -G ${refseq_gtf} -o ${result_path}/${out_name}/${out_name}.gtf ${bam_path}/${out_name}.bam > ${log_path}/${out_name}/${out_name}.Stringtie.log 2>&1
else
    echo "The experiment layout should be either 'single' end or 'pair' end."
fi

# 2.reads distribution -----------------------------------------------------------------------------
echo "2.reads distribution..."

echo "112 `get_time` samtools index ..."
samtools index -@ ${task_num} ${bam_path}/${out_name}.bam

echo "115 `get_time` read_distribution ..."
/usr/bin/python3.8 ${dir_file_current}/read_distribution.py -i ${bam_path}/${out_name}.bam -r ${refseq_bed} > ${log_path}/${out_name}/${out_name}.readsDistribution

