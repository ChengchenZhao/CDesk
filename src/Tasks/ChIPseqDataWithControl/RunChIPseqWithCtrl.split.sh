#!/bin/bash

# Description : ChIP-seq data analysis pipeline.
# Author      : Chengchen Zhao
# Version     : 1.0
# Time        : 2021-07-13 13:29:19

###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################

data_path=${1} # /data/zhaochengchen/bin/pipeline
result_path=${2} # /data/zhaochengchen/bin/pipeline/ChIPseq
out_name=${3} # ChIPseq_test_single

layout=${4} # single
fastq_file=${5} # "YourFastqFile.fastq.gz" for single-end
                # "YourFastqFile_R1.fastq.gz,YourFastqFile_R2.fastq.gz" for pair-end

mapping_index=${6} # /data/zhaochengchen/Data/mm9/mm9
effective_genome_size=${7} # mm
promoter_file=${8} # /data/zhaochengchen/Data/mm9/mm9.refseq.Promoter
chromInfo=${9} # /data/zhaochengchen/Data/mm9/mm9.len

task_num=$cdesk_task_num
if [ ! $cdesk_task_num ]; then
    task_num=8
fi

# data_path=/data/zhaochengchen/Work/OTHER/00.ZCC/H3K9me3/data/GSE89058_hgESCs_H3K9me3
# result_path=/data/zhaochengchen/Work/OTHER/00.ZCC/H3K9me3/data/GSE89058_hgESCs_H3K9me3
# out_name=K9_hiFT
# layout=pair
# fastq_file=GSM2580281_1.fastq,GSM2580281_2.fastq
# mapping_index=/data/zhaochengchen/Data/T2T-CHM13/T2T-CHM13
# effective_genome_size=hg
# promoter_file=/data/zhaochengchen/Data/hg38/hg38.Promoter.bed
# chromInfo=/data/zhaochengchen/Data/hg38/hg38.len
echo "--------------------------------------------------------------------------------" 
if [ ! -n "${9}" ]
then 
    echo "Parameters are missing, please check your input." 
    exit 2
else 
    echo "ChIP-seq data analysis pipeline is now running... ${task_num} ${data_path} ${result_path} ${out_name} ${layout} ${fastq_file} ${mapping_index} ${effective_genome_size} ${promoter_file} ${chromInfo}" 
fi 

out_name=${out_name}.split
log_path=${result_path}/log
bam_path=${result_path}/bam
bw_path=${result_path}/bw
peak_path=${result_path}/peak

file_current=`realpath ${BASH_SOURCE}`
dir_file_current=`dirname ${file_current}`

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
if [ ! -d "${bw_path}" ]; then
  mkdir ${bw_path}
fi
if [ ! -d "${peak_path}" ]; then
  mkdir ${peak_path}
fi
if [ ! -d "${result_path}/${out_name}" ]; then
  mkdir ${result_path}/${out_name}
fi
if [ ! -d "${log_path}/${out_name}" ]; then
  mkdir ${log_path}/${out_name}
fi
if [ ! -d "${peak_path}/${out_name}" ]; then
  mkdir ${peak_path}/${out_name}
fi


get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}
echo ${out_name}
# 1.mapping ----------------------------------------------------------------------------------------
echo "82 `get_time` 1.mapping..."
if [ "${layout}" == "single" ]
then
    echo "85 `get_time` trim_galore ..."
    trim_galore --gzip -j ${task_num} ${data_path}/${fastq_file} --trim-n -o ${data_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    
    echo "88 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -U ${data_path}/${out_name}_trimmed.fq.gz -S ${result_path}/${out_name}/${out_name}.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${data_path}/${out_name}_trimmed.fq.gz

    echo "92 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.sam > ${result_path}/${out_name}/${out_name}.bam
    
    echo "95 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    
    rm -rf ${result_path}/${out_name}/${out_name}.sam ${result_path}/${out_name}/${out_name}.bam

    # Remove duplicate reads
    samtools view -@ ${task_num} -H ${bam_path}/${out_name}.bam > ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -f 0x2 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM") print forward"\n"$0}' >> ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -bS ${bam_path}/${out_name}.Filtered.sam > ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.sam 
    samtools sort -@ ${task_num} -n -o ${bam_path}/${out_name}.Sort.bam ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.bam
    samtools fixmate -@ ${task_num} -m ${bam_path}/${out_name}.Sort.bam ${bam_path}/${out_name}.FixMate.bam
    rm ${bam_path}/${out_name}.Sort.bam
    samtools sort -@ ${task_num} -o ${bam_path}/${out_name}.FixMateSorted.bam ${bam_path}/${out_name}.FixMate.bam
    rm ${bam_path}/${out_name}.FixMate.bam
    samtools markdup -@ ${task_num} -r ${bam_path}/${out_name}.FixMateSorted.bam ${bam_path}/${out_name}.bam
    rm ${bam_path}/${out_name}.FixMateSorted.bam
    
    # Generate bigwig file
    echo "95 `get_time` samtools sort ..."
    macs3 callpeak --SPMR -f BAMPE -t ${bam_path}/${out_name}.bam --bdg -n ${result_path}/${out_name}/${out_name} > ${log_path}/${out_name}/${out_name}.generatebw.log 2>&1
    
    echo "95 `get_time` bedClip ..."
    bedClip ${result_path}/${out_name}/${out_name}_treat_pileup.bdg ${chromInfo} ${bam_path}/${out_name}.tmp.bdg
    
    echo "143 `get_time` bedGraphToBigWig ..."
    grep -v _ ${bam_path}/${out_name}.tmp.bdg | sort -k1,1 -k2,2n  > ${bam_path}/${out_name}.clean.bdg
    bedGraphToBigWig ${bam_path}/${out_name}.clean.bdg ${chromInfo} ${bw_path}/${out_name}.bw
    rm -rf ${result_path}/${out_name}/${out_name}_treat_pileup.bdg
    rm -rf ${result_path}/${out_name}/${out_name}_control_lambda.bdg
    rm -rf ${result_path}/${out_name}/${out_name}_peaks.narrowPeak
    rm -rf ${result_path}/${out_name}/${out_name}_peaks.xls
    rm -rf ${result_path}/${out_name}/${out_name}_summits.bed
    rm -rf ${bam_path}/${out_name}.tmp.bdg
    rm -rf ${bam_path}/${out_name}.clean.bdg
elif [ "${layout}" == "pair" ]
then
    echo "133 `get_time` trim_galore ..."
    oldIFS=$IFS
    IFS=','
    read -ra fastq_file_arr1 <<< "${fastq_file}"
    trim_galore --gzip -j ${task_num} --paired ${data_path}/${fastq_file_arr1[0]} ${data_path}/${fastq_file_arr1[1]} --trim-n -o ${data_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    IFS=$oldIFS

    echo "138 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -1 ${data_path}/${out_name}_val_1.fq.gz -2 ${data_path}/${out_name}_val_2.fq.gz -S ${result_path}/${out_name}/${out_name}.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${data_path}/${out_name}_val_1.fq.gz ${data_path}/${out_name}_val_2.fq.gz
    
    echo "142 `get_time` samtools ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.sam > ${result_path}/${out_name}/${out_name}.bam
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    rm -rf ${result_path}/${out_name}/${out_name}.sam ${result_path}/${out_name}/${out_name}.bam

    # Remove duplicate reads
    samtools view -@ ${task_num} -H ${bam_path}/${out_name}.bam > ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -f 0x2 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM") print forward"\n"$0}' >> ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -bS ${bam_path}/${out_name}.Filtered.sam > ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.sam 
    samtools sort -@ ${task_num} -n -o ${bam_path}/${out_name}.Sort.bam ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.bam
    samtools fixmate -@ ${task_num} -m ${bam_path}/${out_name}.Sort.bam ${bam_path}/${out_name}.FixMate.bam
    rm ${bam_path}/${out_name}.Sort.bam
    samtools sort -@ ${task_num} -o ${bam_path}/${out_name}.FixMateSorted.bam ${bam_path}/${out_name}.FixMate.bam
    rm ${bam_path}/${out_name}.FixMate.bam
    samtools markdup -@ ${task_num} -r ${bam_path}/${out_name}.FixMateSorted.bam ${bam_path}/${out_name}.bam
    rm ${bam_path}/${out_name}.FixMateSorted.bam
    
    # Generate bigwig file
    echo "162 `get_time` macs3 callpeak ..."
    macs3 callpeak --SPMR -f BAMPE -t ${bam_path}/${out_name}.bam --bdg -n ${result_path}/${out_name}/${out_name} > ${log_path}/${out_name}/${out_name}.generatebw.log 2>&1

    # bug fix: bw generate error
    grep -v "-" ${result_path}/${out_name}/${out_name}_treat_pileup.bdg > ${bam_path}/${out_name}.tmp.bdg
    mv ${bam_path}/${out_name}.tmp.bdg ${result_path}/${out_name}/${out_name}_treat_pileup.bdg

    echo "165 `get_time` bedClip ..."
    bedClip ${result_path}/${out_name}/${out_name}_treat_pileup.bdg ${chromInfo} ${bam_path}/${out_name}.tmp.bdg

    echo "168 `get_time` bedGraphToBigWig ..."
    grep -v _ ${bam_path}/${out_name}.tmp.bdg | sort -k1,1 -k2,2n > ${bam_path}/${out_name}.clean.bdg
    bedGraphToBigWig ${bam_path}/${out_name}.clean.bdg ${chromInfo} ${bw_path}/${out_name}.bw
    rm -rf ${result_path}/${out_name}/${out_name}_treat_pileup.bdg
    rm -rf ${result_path}/${out_name}/${out_name}_control_lambda.bdg
    rm -rf ${result_path}/${out_name}/${out_name}_peaks.narrowPeak
    rm -rf ${result_path}/${out_name}/${out_name}_peaks.xls
    rm -rf ${result_path}/${out_name}/${out_name}_summits.bed
    rm -rf ${bam_path}/${out_name}.tmp.bdg
    rm -rf ${bam_path}/${out_name}.clean.bdg

else
    echo "The experiment layout should be either 'single' end or 'pair' end."
fi