#!/bin/bash

# Description : Preprocessing+Mapping+QC
# Author      : WEI SHI
# Version     : 1.0
# Time        : 2024-10-19 22:40:00

show_help() {
    echo "Usage: BulkRNAseqData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species (mm9, mm10)"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (fastq.gz), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -l 1 or 2		1:Single sequencing, 2:Pair sequencing (default is 2)"
}

THREAD=8
LAYOUT=2

while getopts ":hs:i:o:p:l:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            INPUT_DIRECTORY=${OPTARG}
            ;;
        s )
            SPECIES=${OPTARG}
            ;;
        o)
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        p )
            THREAD=${OPTARG}
            ;;
        l )
            LAYOUT=${OPTARG}
            if [[ "$LAYOUT" != "1" && "$LAYOUT" != "2" ]]; then
                echo "Invalid value for -l: $LAYOUT. It must be 1 or 2." 1>&2
                show_help
                exit 1
            fi
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            show_help
            exit 1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            show_help
            exit 1
            ;;
    esac
done

# 检查是否提供了必要的参数
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$SPECIES" ] || [ -z "$OUTPUT_DIRECTORY" ]; then
    echo "Error: INPUT_DIRECTORY, OUTPUT_DIRECTORY and SPECIES must be provided."
    show_help
    exit 1
fi

# 分配参考数据
if [ "$SPECIES" = "mm9" ]; then
    mapping_index="/mnt/shiwei/Ref/mm9/mm9"
    mapping_gtf="/mnt/shiwei/Ref/mm9/mm9_refseq.gtf"
    gfold_gtf="/mnt/shiwei/Ref/mm9/mm9.ncbiRefSeq.gtf"
    RSeQC_bed="/mnt/shiwei/Ref/mm9/mm9.refseq.bed"
elif [ "$SPECIES" = "mm10" ]; then
    mapping_index="/mnt/shiwei/Ref/mm10/mm10"
    mapping_gtf="/mnt/shiwei/Ref/mm10/mm10_refseq.gtf"
    gfold_gtf="/mnt/shiwei/Ref/mm10/mm10.ncbiRefSeq.gtf"
    RSeQC_bed="/mnt/shiwei/Ref/mm10/mm10.refseq.bed"
fi

# 提取文件名
if [ "$LAYOUT" = "1" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.fastq.gz | sed 's/.fastq.gz//' | sort -u)
elif [ "$LAYOUT" = "2" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.fastq.gz | sed 's/_1.fastq.gz//;s/_2.fastq.gz//' | sort -u)
fi


log_path="${OUTPUT_DIRECTORY}/Log"
bam_path="${OUTPUT_DIRECTORY}/Bam"
QC_path="${OUTPUT_DIRECTORY}/QC"
Expression_path="${OUTPUT_DIRECTORY}/Expression"


if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir -p ${OUTPUT_DIRECTORY}
fi
if [ ! -d "${log_path}" ]; then
    mkdir -p ${log_path}
fi
if [ ! -d "${bam_path}" ]; then
    mkdir -p ${bam_path}
fi
if [ ! -d "${QC_path}" ]; then
    mkdir -p ${QC_path}
fi
if [ ! -d "${Expression_path}" ]; then
    mkdir -p ${Expression_path}
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}


####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################
counter=0

for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    echo "-----------------------------------------------------------------------------------------------------"
    echo "----------------------------第 ${counter} 个样本，样本名为 ${SAMPLE_PREFIX}--------------------------"
    echo "-----------------------------------------------------------------------------------------------------"

    FILE_NAME=""
    R1_FILE_NMAE=""
    R2_FILE_NMAE=""

    if [ "${LAYOUT}" == "1" ]; then
        FILE_NAME="${FILE}.fastq.gz"
    elif [ "${LAYOUT}" == "2" ]; then
        R1_FILE_NMAE="${FILE}_1.fastq.gz"
        R2_FILE_NMAE="${FILE}_2.fastq.gz"
    fi

    echo " "
    echo "--------------------------------------------INITIALIZING----------------------------------------------" 
    echo "RNA-seq data analysis pipeline is now running..."
    echo "Number of threads ---------- ${THREAD}" 
    echo "Directory of data ---------- ${INPUT_DIERCTORY}"
    echo "Directory of result ---------- ${OUTPUT_DIERCTORY}"
    echo "Name of sample ---------- ${SAMPLE_PREFIX}"
    echo "File of data ---------- ${FILE_NAME}"
    echo "File of R1 data ---------- ${R1_FILE_NMAE}"
    echo "File of R2 data ---------- ${R2_FILE_NMAE}"
    echo "Method of sequencing ---------- ${LAYOUT}"
    echo "Mapping index ---------- ${mapping_index}"
    echo "Mapping gtf ---------- ${mapping_gtf}"
    echo "GFOLD gtf ---------- ${gfold_gtf}"
    echo "RSeQC bed ---------- ${RSeQC_bed}"
    echo "----------------------------------INITIALIZATION FINISHED-------------------------------------"
    echo " "

    echo " "
    echo "-----------------------------------GENERATING DIRECTORY-------------------------------------" 
    echo "Please wait, the directory for storing the results is being generated..."


    echo "Save LOG records for all samples ---------- ${log_path}" 
    echo "Store BAM files for all samples ---------- ${bam_path}"
    echo "Store quality control results, including fastqc results, mapping results, and RSeQC results ---------- ${QC_path}"
    echo "Store gene expression data, including count and fpkm versions ---------- ${Expression_path}"
    echo "---------------------------DIRECTORY GENERATION FINISHED-----------------------------"
    echo " "

    echo " "
    echo "---------------------------------------START ANALYZING-----------------------------------------" 
    echo "Preparation work completed, start analysis..."

    # 1.mapping ----------------------------------------------------------------------------------------
    echo " "
    echo "1.mapping..."
    echo " "
    if [ "${LAYOUT}" == "1" ]; then
        echo "67 `get_time` trim_galore ..."
        echo " "
	trim_galore --gzip -j ${THREAD} ${FILE_NAME} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
    elif [ "${LAYOUT}" == "2" ]; then
        echo "86 `get_time` trim_galore ..."
        echo " "
	trim_galore --gzip -j ${THREAD} --paired ${R1_FILE_NMAE} ${R2_FILE_NMAE} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
    fi


    if [ "${LAYOUT}" == "1" ]; then
        echo "70 `get_time` hisat2 ..."
        echo " "
	hisat2 -p ${THREAD} --dta -x ${mapping_index} -U ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_trimmed.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1 
        rm -rf ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_trimmed.fq.gz
    elif [ "${LAYOUT}" == "2" ]; then
        echo "91 `get_time` hisat2 ..."
        echo " "
	hisat2 -p ${THREAD} --dta -x ${mapping_index} -1 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_1.fq.gz -2 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_2.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
        rm -rf ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_1.fq.gz ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_2.fq.gz
    fi



    echo "74 & 95 `get_time` samtools view ..."
    echo " "
    samtools view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${bam_path}/${SAMPLE_PREFIX}.bam


    echo "77 & 98 `get_time` samtools sort ..."
    echo " "
    samtools sort -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam > ${bam_path}/${SAMPLE_PREFIX}.Sorted.bam
    rm -rf ${bam_path}/${SAMPLE_PREFIX}.align.sam ${bam_path}/${SAMPLE_PREFIX}.bam
    mv ${bam_path}/${SAMPLE_PREFIX}.Sorted.bam ${bam_path}/${SAMPLE_PREFIX}.bam

    echo "112 `get_time` samtools index ..."
    echo " "
    samtools index -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam



    # 2.Calculate gene expression levels --------------------------------------------------------------------
    echo " "
    echo "2.Calculating gene expression levels ..."
    echo " "
    echo "116 `get_time` gfold count ..."
    echo " "
    samtools view -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam | gfold count -ann ${gfold_gtf} -tag stdin -o ${Expression_path}/${SAMPLE_PREFIX}.read_cnt


    if [ ! -d "${Expression_path}/Stringtie" ]; then
        mkdir -p ${Expression_path}/Stringtie
    fi

    if [ ! -d "${bam_path}/Stringtie" ]; then
        mkdir -p ${bam_path}/Stringtie
    fi

    
    echo "160 `get_time` stringtie count ..."
    echo " "
    nohup stringtie -e -B -p ${THREAD} -G ${mapping_gtf} -o ${Expression_path}/Stringtie/${SAMPLE_PREFIX}.gtf ${bam_path}/${SAMPLE_PREFIX}.bam > ${log_path}/${SAMPLE_PREFIX}.Stringtie.log 2>&1


    # 3.QC --------------------------------------------------------------------
    echo " "
    echo "3.QC ..."
    echo " "
    echo "120 `get_time` fastqc ..."
    echo " "
    fastqc -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz

    echo "124 `get_time` RSeQC ..."
    echo " "
    read_distribution.py -i ${bam_path}/${SAMPLE_PREFIX}.bam -r ${RSeQC_bed} > ${QC_path}/${SAMPLE_PREFIX}_dirstribution.txt

    echo "---------------------------ANALYZING FINISHED-----------------------------"
    echo " "

done

echo " "
echo "begin merge fpkm expression..."
echo " "
# -------------------------------------------合并fpkm----------------------------------------------
# 提取第一个.read_cnt文件的第一列并追加到CSV文件中
awk '{print $1}' $(ls ${Expression_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# 遍历所有的.read_cnt文件
for file in ${Expression_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $5}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

echo "gene_name,$(ls ${Expression_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv

# 将原始数据追加到临时文件中
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv
# 用临时文件替换原始文件
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_fpkm.csv

rm -rf ${Expression_path}/output1.csv
echo " "
echo "fpkm expression has been merged..."
echo " "

echo " "
echo "begin merge counts expression..."
echo " "
# -------------------------------------------合并count----------------------------------------------
# 提取第一个.read_cnt文件的第一列并追加到CSV文件中
awk '{print $1}' $(ls ${Expression_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# 遍历所有的.read_cnt文件
for file in ${Expression_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $3}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

echo "gene_name,$(ls ${Expression_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv

# 将原始数据追加到临时文件中
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv
# 用临时文件替换原始文件
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_count.csv

rm -rf ${Expression_path}/output1.csv
echo " "
echo "counts expression has been merged..."
echo " "




