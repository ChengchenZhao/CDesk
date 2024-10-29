#!/bin/bash

# Description : Calculate the similarity between the collected sample and a set of standard data
# Author      : WEI SHI
# Version     : 1.0
# Time        : 2024-10-22 13:23:00

show_help() {
    echo "Usage: bash BulkRNAseqSim.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help          	Show this help message"
    echo "  -b BASE_INPUT_FILE       	Specify the base_input file (fpkm.csv), column name as sample, row name as gene"
    echo "  -x BASE_REPLICATE       	Specify the symbol representing repeated experiments of base data (if not, fill in NO)"
    echo "  -r REF_INPUT_FILE       	Specify the ref_input file (fpkm.csv), column name as sample, row name as gene"
    echo "  -y REF_REPLICATE       	Specify the symbol representing repeated experiments of reference data (if not, fill in NO)"
    echo "  -o OUTPUT_DIRECTORY 	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -g GENE_LIST  		Specify the gene list (if not, fill in ALL)"
}


while getopts ":hb:x:r:y:o:g:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        b )
            BASE_INPUT_FILE=${OPTARG}
            ;;
        x )
            BASE_REPLICATE=${OPTARG}
            ;;
        r )
            REF_INPUT_FILE=${OPTARG}
            ;;
        y )
            REF_REPLICATE=${OPTARG}
            ;;
        o )
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        g )
            GENE_LIST=${OPTARG}
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
if [ -z "$BASE_INPUT_FILE" ] || [ -z "$REF_INPUT_FILE" ] || [ -z "$OUTPUT_DIRECTORY" ] || [-z "$BASE_REPLICATE"] || [-z "$REF_REPLICATE"] || [-z "$GENE_LIST"]; then
    echo "Error: BASE_INPUT_FILE, BASE_REPLICATE, REF_INPUT_FILE, REF_REPLICATE, OUTPUT_DIRECTORY and GENE_LIST must be provided."
    show_help
    exit 1
fi

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir ${OUTPUT_DIRECTORY}
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

echo "BASE_INPUT_FILE:		${BASE_INPUT_FILE}"
echo "BASE_REPLICATE:          ${BASE_REPLICATE}"
echo "REF_INPUT_FILE:          ${REF_INPUT_FILE}"
echo "REF_REPLICATE:          ${REF_REPLICATE}"
echo "OUTPUT_DIRECTORY:          ${OUTPUT_DIRECTORY}"
echo "GENE_LIST:          ${GENE_LIST}"


Rscript /mnt/shiwei/bin/Similarity.R ${BASE_INPUT_FILE} ${BASE_REPLICATE} ${REF_INPUT_FILE} ${REF_REPLICATE} ${OUTPUT_DIRECTORY} ${GENE_LIST}
echo "136 `get_time` Processing completed ..."


