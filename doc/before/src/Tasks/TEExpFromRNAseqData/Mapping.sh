#!/bin/bash

# Description : ATAC-seq data analysis pipeline.
# Author      : Chengchen Zhao
# Version     : 1.1
# Time        : 2021-07-27 13:26:10

####################################################################################################
#################################     PARAMETER/FUNCTION DEFINE     ################################
####################################################################################################
data_path=${1} # /data/zhaochengchen/bin/pipeline
result_path=${2} # /data/zhaochengchen/bin/pipeline/ATAC
out_name=${3} # ATAC_test_single

layout=${4} # single/pair
fastq_file=${5} # "YourFastqFile.fastq.gz" for single-end
                # "YourFastqFile_R1.fastq.gz,YourFastqFile_R2.fastq.gz" for pair-end
mapping_index=${6} # /data/zhaochengchen/Data/mm9/mm9
species=${7} # mm
promoter_file=${8} # /data/zhaochengchen/Data/mm9/mm9.refseq.Promoter
chromInfo=${9} # /data/zhaochengchen/Data/mm9/mm9.len

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
    echo "ATAC-seq data analysis pipeline is now running... ${task_num} ${data_path}  ${result_path}  ${out_name}  ${layout}  ${fastq_file}  ${mapping_index}  ${species}  ${promoter_file}  ${chromInfo}" 
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
    # echo "74 `get_time` trim_galore ..."
    # trim_galore --gzip -j ${task_num} ${data_path}/${fastq_file} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    
    echo "77 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -U ${data_path}/${fastq_file} -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    # rm -rf ${result_path}/${out_name}_trimmed.fq.gz

    echo "81 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam

    echo "84 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam

elif [ "${layout}" == "pair" ]
then
    oldIFS=$IFS
    IFS=','
    read -ra fastq_file_arr <<< "${fastq_file}"
    echo "132 `get_time` trim_galore ..."
    # trim_galore --gzip -j ${task_num} --paired ${data_path}/${fastq_file_arr[0]} ${data_path}/${fastq_file_arr[1]} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    

    echo "135 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -1 ${data_path}/${fastq_file_arr[0]} -2 ${data_path}/${fastq_file_arr[1]} -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    # rm -rf ${result_path}/${out_name}_val_1.fq.gz ${result_path}/${out_name}_val_2.fq.gz
    
    IFS=$oldIFS

    echo "139 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam

    echo "142 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam

else
    echo "The experiment layout should be either 'single' end or 'pair' end."
fi

echo "Done!"
echo "--------------------------------------------------------------------------------" 
