#!/bin/bash

# Description : ChIP-seq data analysis pipeline.
# Author      : Chengchen Zhao
# Version     : 1.0
# Time        : 2021-07-13 13:29:19

###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################

result_path=${1} # /data/zhaochengchen/bin/pipeline/ChIPseq
out_name=${2} # ChIPseq_test_single
ctrl_out_name=${3} # ChIPseq_test_single

effective_genome_size=${4} # mm
promoter_file=${5} # /data/zhaochengchen/Data/mm9/mm9.refseq.Promoter
chromInfo=${6} # /data/zhaochengchen/Data/mm9/mm9.len

shiftsize=${7}


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
if [ ! -n "${7}" ]
then 
    echo "Parameters are missing, please check your input." 
    exit 2
else 
    echo "ChIP-seq data analysis pipeline is now running... ${task_num} ${result_path} ${out_name} ${layout} ${ctrl_out_name} ${effective_genome_size} ${promoter_file} ${chromInfo} ${treat_fastq_file}" 
fi 

log_path=${result_path}/log
bam_path=${result_path}/bam
bw_path=${result_path}/bw
peak_path=${result_path}/peak

file_current=`realpath ${BASH_SOURCE}`
dir_file_current=`dirname ${file_current}`

####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################

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


# 1.Generate bigwig file -----------------------------------------------------------------------------
echo "155 `get_time` BigWig ..."
TreatReadsNumber=`samtools view -c -F 260 ${bam_path}/${out_name}.split.bam`
CtrlReadsNumber=`samtools view -c -F 260 ${bam_path}/${ctrl_out_name}.split.bam`
scaleFactor1=`bc -l <<< "1000000/${TreatReadsNumber}"`
scaleFactor2=`bc -l <<< "1000000/${CtrlReadsNumber}"`
bigwigCompare -b1 ${bw_path}/${out_name}.split.bw -b2 ${bw_path}/${ctrl_out_name}.split.bw --scaleFactors ${scaleFactor1}:${scaleFactor2} --skipZeroOverZero --operation log2 -bs 10 -o ${bw_path}/${out_name}.bw -of bigwig
bigWigToBedGraph ${bw_path}/${out_name}.bw ${bw_path}/${out_name}.tmp.bdg
awk -v "OFS=\t" '{if ($4<0) print ($1,$2,$3,0);else print ($0)}' ${bw_path}/${out_name}.tmp.bdg > ${bw_path}/${out_name}.bdg
bedGraphToBigWig ${bw_path}/${out_name}.bdg ${chromInfo} ${bw_path}/${out_name}.bw
rm ${bw_path}/${out_name}.tmp.bdg ${bw_path}/${out_name}.bdg

# 2.peak calling -----------------------------------------------------------------------------
echo "184 `get_time` 2.peak calling..."
echo "185 `get_time` macs14 callpeak ..."
macs14 callpeak -g ${effective_genome_size} -p 1e-3 --keep-dup 1 --nomodel --shiftsize ${shiftsize} -t ${bam_path}/${out_name}.split.bam -c ${bam_path}/${ctrl_out_name}.split.bam -n ${peak_path}/${out_name}/${out_name}.e3 > ${result_path}/log/${out_name}/${out_name}.callpeak1e_3.log 2>&1
macs14 callpeak -g ${effective_genome_size} -p 1e-5 --keep-dup 1 --nomodel --shiftsize ${shiftsize} -t ${bam_path}/${out_name}.split.bam -c ${bam_path}/${ctrl_out_name}.split.bam -n ${peak_path}/${out_name}/${out_name}.e5 > ${result_path}/log/${out_name}/${out_name}.callpeak1e_5.log 2>&1
macs14 callpeak -g ${effective_genome_size} -p 1e-7 --keep-dup 1 --nomodel --shiftsize ${shiftsize} -t ${bam_path}/${out_name}.split.bam -c ${bam_path}/${ctrl_out_name}.split.bam -n ${peak_path}/${out_name}/${out_name}.e7 > ${result_path}/log/${out_name}/${out_name}.callpeak1e_7.log 2>&1

# 3.Promoter signal calculation ----------------------------------------------------------------------
echo "189 `get_time` 3.Promoter signal calculation..."
echo "190 `get_time` GetBigwigSignal ..."
python3 ${dir_file_current}/GetBigwigSignal.py ${promoter_file} ${bw_path}/${out_name}.bw ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp 1 > ${log_path}/${out_name}/${out_name}.GetBigwigSignal.log 2>&1
cut -f 5 ${promoter_file} | paste - ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp > ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2

echo "194 `get_time` GeneSymbolUniq ..."
python3 ${dir_file_current}/GeneSymbolUniq.py ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2 ${result_path}/${out_name}/${out_name}.PromoterSignal 0
rm -rf ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp
rm -rf ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2

echo "199 `get_time` 4.Summary..."
echo "200 `get_time` ChIPseqMappingSummary ..."
python ${dir_file_current}/ChIPseqMappingSummary.py ${result_path}/log/${out_name}/${out_name}.Mapping.log ${result_path}/MappingSummary_ChIPseq.txt ${out_name}
wc -l ${peak_path}/${out_name}/${out_name}.e3_peaks.bed >> ${result_path}/PeakNumber_1e3_ChIPseq.txt
wc -l ${peak_path}/${out_name}/${out_name}.e5_peaks.bed >> ${result_path}/PeakNumber_1e5_ChIPseq.txt
wc -l ${peak_path}/${out_name}/${out_name}.e7_peaks.bed >> ${result_path}/PeakNumber_1e7_ChIPseq.txt

echo "207 `get_time` Done!"
echo "--------------------------------------------------------------------------------" 
