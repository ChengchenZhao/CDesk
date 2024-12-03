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
effective_genome_size=${7} # mm
promoter_file=${8} # /data/zhaochengchen/Data/mm9/mm9.refseq.Promoter
chromInfo=${9} # /data/zhaochengchen/Data/mm9/mm9.len
shiftsize=${10}

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
    echo "ATAC-seq data analysis pipeline is now running... ${task_num} ${data_path}  ${result_path}  ${out_name}  ${layout}  ${fastq_file}  ${mapping_index}  ${effective_genome_size}  ${promoter_file}  ${chromInfo}" 
fi 

log_path=${result_path}/log
bam_path=${result_path}/bam
bw_path=${result_path}/bw
peak_path=${result_path}/peak

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
if [ ! -d "${result_path}/FragmentsLength" ]; then
  mkdir ${result_path}/FragmentsLength
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

# 1.mapping ----------------------------------------------------------------------------------------
echo "1.mapping..."
if [ "${layout}" == "single" ]
then
    echo "74 `get_time` trim_galore ..."
    trim_galore --gzip -j ${task_num} ${data_path}/${fastq_file} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    
    echo "77 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -U ${result_path}/${out_name}_trimmed.fq.gz -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${result_path}/${out_name}_trimmed.fq.gz

    echo "81 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam

    echo "84 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam

    echo "88 `get_time` chromosome distribution ..."
    # chromosome distribution
    samtools view -@ ${task_num} ${bam_path}/${out_name}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${result_path}/${out_name}/${out_name}.ChromosomeDistribution.txt
    samtools view -@ ${task_num} -F 0x400 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr") print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${result_path}/${out_name}/${out_name}.Filtered.ChromosomeDistribution.txt
    

    echo "94 `get_time` filtering ..."
    # filtering
    samtools view -@ ${task_num} -H ${bam_path}/${out_name}.bam > ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -F 0x400 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM") print forward"\n"$0}' >> ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -bS ${bam_path}/${out_name}.Filtered.sam > ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.sam 

    echo "101 `get_time` bamToBed ..."
    bamToBed -i ${bam_path}/${out_name}.Filtered.bam > ${result_path}/${out_name}/${out_name}.UniqSEReads.bed
    n=`wc -l ${result_path}/${out_name}/${out_name}.UniqSEReads.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"` 

    echo "106 `get_time` genomeCoverageBed ..."
    genomeCoverageBed -bga -scale $c -i ${result_path}/${out_name}/${out_name}.UniqSEReads.bed -g ${chromInfo} > ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg 

    echo "109 `get_time` bedClip ..."
    bedClip ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg ${chromInfo} ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg

    echo "112 `get_time` bedSort ..."
    bedSort ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp

    echo "115 `get_time` bedGraphToBigWig ..."
    bedGraphToBigWig ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp ${chromInfo} ${bw_path}/${out_name}.bw
    rm ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg
    rm ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp
    rm ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg

    echo "121 `get_time` bedSort ..."
    bedSort ${result_path}/${out_name}/${out_name}.UniqSEReads.bed ${result_path}/${out_name}/${out_name}.UniqSEReadsSorted.bed

    echo "124 `get_time` bedToBigBed ..."
    bedToBigBed ${result_path}/${out_name}/${out_name}.UniqSEReadsSorted.bed ${chromInfo} ${result_path}/${out_name}/${out_name}.UniqSEReads.bd

    ln -s ${result_path}/${out_name}/${out_name}.UniqSEReads.bed ${result_path}/${out_name}/${out_name}.Fragments.bed

elif [ "${layout}" == "pair" ]
then
    oldIFS=$IFS
    IFS=','
    read -ra fastq_file_arr <<< "${fastq_file}"
    echo "132 `get_time` trim_galore ..."
    trim_galore --gzip -j ${task_num} --paired ${data_path}/${fastq_file_arr[0]} ${data_path}/${fastq_file_arr[1]} --trim-n -o ${result_path} --no_report_file --basename ${out_name} > ${log_path}/${out_name}/${out_name}.TrimGalore.log 2>&1
    IFS=$oldIFS
    echo "135 `get_time` bowtie2 ..."
    bowtie2 -q -p ${task_num} -x ${mapping_index} -1 ${result_path}/${out_name}_val_1.fq.gz -2 ${result_path}/${out_name}_val_2.fq.gz -S ${result_path}/${out_name}/${out_name}.align.sam > ${log_path}/${out_name}/${out_name}.Mapping.log 2>&1
    rm -rf ${result_path}/${out_name}_val_1.fq.gz ${result_path}/${out_name}_val_2.fq.gz
    
    echo "139 `get_time` samtools view ..."
    samtools view -@ ${task_num} -bS ${result_path}/${out_name}/${out_name}.align.sam > ${result_path}/${out_name}/${out_name}.bam

    echo "142 `get_time` samtools sort ..."
    samtools sort -@ ${task_num} ${result_path}/${out_name}/${out_name}.bam > ${bam_path}/${out_name}.bam
    rm -rf ${result_path}/${out_name}/${out_name}.align.sam ${result_path}/${out_name}/${out_name}.bam

    # chromosome distribution
    echo "147 `get_time` samtools view ..."
    samtools view -@ ${task_num} ${bam_path}/${out_name}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${result_path}/${out_name}/${out_name}.ChromosomeDistribution.txt
    samtools view -@ ${task_num} -f 0x2 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr") print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${result_path}/${out_name}/${out_name}.Filtered.ChromosomeDistribution.txt
    
    # filtering
    echo "152 `get_time` samtools view ..."
    samtools view -@ ${task_num} -H ${bam_path}/${out_name}.bam > ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -f 0x2 ${bam_path}/${out_name}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM") print forward"\n"$0}' >> ${bam_path}/${out_name}.Filtered.sam
    samtools view -@ ${task_num} -bS ${bam_path}/${out_name}.Filtered.sam > ${bam_path}/${out_name}.Filtered.bam
    rm ${bam_path}/${out_name}.Filtered.sam 

    echo "158 `get_time` bamToBed ..."
    bamToBed -i ${bam_path}/${out_name}.Filtered.bam > ${result_path}/${out_name}/${out_name}.RawReads.bed

    echo "161 `get_time` samtools sort ..."
    samtools sort -n -@ ${task_num} ${bam_path}/${out_name}.Filtered.bam > ${bam_path}/${out_name}.FilteredSorted.bam
    rm ${bam_path}/${out_name}.Filtered.bam

    echo "165 `get_time` bamToBed ..."
    bamToBed -bedpe -i ${bam_path}/${out_name}.FilteredSorted.bam | awk '{print $1"\t"$2"\t"$6}' | sort -k1,1 -k2,2n | uniq | grep -v WARNING > ${result_path}/${out_name}/${out_name}.Fragments.bed 2>&1
    awk '{print $3-$2}' ${result_path}/${out_name}/${out_name}.Fragments.bed | sort | uniq -c | sort -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > ${result_path}/${out_name}/${out_name}.FragmentsLength.txt
    sort -k1,1 -k2,2g ${result_path}/${out_name}/${out_name}.Fragments.bed > ${result_path}/${out_name}/${out_name}.UniqSEReads.bed
    
    echo "170 `get_time` grep chrEBV ..."
    cat ${result_path}/${out_name}/${out_name}.UniqSEReads.bed | grep -v chrEBV > ${result_path}/${out_name}/${out_name}.UniqSEReadsGrep.bed
    mv ${result_path}/${out_name}/${out_name}.UniqSEReadsGrep.bed ${result_path}/${out_name}/${out_name}.UniqSEReads.bed
    
    n=`wc -l ${result_path}/${out_name}/${out_name}.UniqSEReads.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"` 
    
    echo "177 `get_time` genomeCoverageBed ..."
    genomeCoverageBed -bga -scale $c -i ${result_path}/${out_name}/${out_name}.UniqSEReads.bed -g ${chromInfo} > ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg
    
    echo "180 `get_time` bedClip ..."
    bedClip ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg ${chromInfo} ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg
    
    echo "183 `get_time` bedSort ..."
    bedSort ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp
    
    echo "186 `get_time` bedGraphToBigWig ..."
    bedGraphToBigWig ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp ${chromInfo} ${bw_path}/${out_name}.bw
    rm ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg
    rm ${result_path}/${out_name}/${out_name}.UniqSEReads.bdg.temp
    rm ${result_path}/${out_name}/${out_name}.UniqSEReadsClip.bdg
    
    echo "192 `get_time` bedSort ..."
    bedSort ${result_path}/${out_name}/${out_name}.UniqSEReads.bed ${result_path}/${out_name}/${out_name}.UniqSEReadsSorted.bed
    
    echo "195 `get_time` bedToBigBed ..."
    bedToBigBed ${result_path}/${out_name}/${out_name}.UniqSEReadsSorted.bed ${chromInfo} ${result_path}/${out_name}/${out_name}.UniqSEReads.bd

else
    echo "The experiment layout should be either 'single' end or 'pair' end."
fi

# 2.peak calling -----------------------------------------------------------------------------
echo "2.peak calling..."
# keep all is beacause the each bed file have been uniqued then merge. The "keep-all" parameter means the reads from different replicates are keeped.
# single-end input : reads in + and reads in - would be shifted by macs14
echo "206 `get_time` macs14 ..."
macs14 -t ${result_path}/${out_name}/${out_name}.Fragments.bed -n ${peak_path}/${out_name}.e3 -f BED -g ${effective_genome_size} --keep-dup all --nomodel --shiftsize ${shiftsize} -p 1e-3 > ${log_path}/${out_name}/${out_name}.e3.PeakCalling.log 2>&1
macs14 -t ${result_path}/${out_name}/${out_name}.Fragments.bed -n ${peak_path}/${out_name}.e5 -f BED -g ${effective_genome_size} --keep-dup all --nomodel --shiftsize ${shiftsize} -p 1e-5 > ${log_path}/${out_name}/${out_name}.e5.PeakCalling.log 2>&1
macs14 -t ${result_path}/${out_name}/${out_name}.Fragments.bed -n ${peak_path}/${out_name}.e7 -f BED -g ${effective_genome_size} --keep-dup all --nomodel --shiftsize ${shiftsize} -p 1e-7 > ${log_path}/${out_name}/${out_name}.e7.PeakCalling.log 2>&1
# cut -f 8 ${peak_path}/${out_name}_peaks.xls | grep -v "#" > ${peak_path}/${out_name}.PeakFoldEnrichment.txt
# Rscript PeakFoldEnrichmentDistribution.r ${peak_path}/${out_name}.PeakFoldEnrichment.txt ${peak_path}/${out_name}
# awk -v "OFS=\t" '{if ($8 > 5 && substr($1,1,3) == "chr" && $1 != "chr") printf ("%s\t%d\t%d\tmacspeak_%d\n",$1,$2-1,$3,NR)}' ${peak_path}/${out_name}_peaks.xls > ${peak_path}/${out_name}.FE5Peaks.bed
# awk -v "OFS=\t" '{if ($8 > 5 && substr($1,1,3) == "chr" && $1 != "chr") printf ("%s\t%d\t%d\tmacspeak_%d\n",$1,$2-2+$5,$2-1+$5,NR)}' ${peak_path}/${out_name}_peaks.xls > ${peak_path}/${out_name}.FE5Summits.bed

# 3.Promoter signal calculation ----------------------------------------------------------------------
echo "3.Promoter signal calculation..."

echo "218 `get_time` GetBigwigSignal.py ..."
/usr/bin/python3.8 ${dir_file_current}/GetBigwigSignal.py ${promoter_file} ${bw_path}/${out_name}.bw ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp 1

echo "221 `get_time` cut ..."
cut -f 5 ${promoter_file} | paste - ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp > ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2

echo "224 `get_time` GeneSymbolUniq.py ..."
/usr/bin/python3.8 ${dir_file_current}/GeneSymbolUniq.py ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2 ${result_path}/${out_name}/${out_name}.PromoterSignal 0
rm ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp
rm ${result_path}/${out_name}/${out_name}.PromoterSignal.tmp2


# 4.Summary ----------------------------------------------------------------------
echo "4.Summary..."

echo "233 `get_time` ATACseqMappingSummary.py ..."
/usr/bin/python3.8 ${dir_file_current}/ATACseqMappingSummary.py ${result_path}/log/${out_name}/${out_name}.Mapping.log ${result_path}/${out_name}/${out_name}.ChromosomeDistribution.txt ${result_path}/MappingSummary_ATACseq.txt ${out_name}

echo "236 `get_time` macs14.py ..."
macs14 -t ${result_path}/${out_name}/${out_name}.UniqSEReadsSorted.bed -n ${peak_path}/${out_name}.HC -f BED -g ${effective_genome_size} --keep-dup all --nomodel --shiftsize ${shiftsize} -p 1e-7 > ${result_path}/log/${out_name}/${out_name}.HighConfidentPeakCalling.log 2>&1

echo "239 `get_time` FragmentLengthDistribution.r ..."
/usr/bin/Rscript ${dir_file_current}/FragmentLengthDistribution.R ${result_path}/${out_name}/${out_name}.FragmentsLength.txt ${out_name} ${result_path}/FragmentsLength & > /dev/null

echo "Done!"
echo "--------------------------------------------------------------------------------" 
