#!/bin/bash

# Description : Test the experimental repeatability of each sample
# Author      : WEI SHI
# Version     : 1.0
# Time        : 2024-10-20 15:29:00

show_help() {
    echo "Usage: bash BulkRNAseqRep.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -i INPUT_FILE	Specify the input file (.csv), column name as sample, row name as gene"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -r REPLICATE_INDEX	Specify the symbol representing repeated experiments."
}

while getopts ":hi:o:r:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            INPUT_DIRECTORY=${OPTARG}
            ;;
        o)
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        r )
            REPLICATE_INDEX=${OPTARG}
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
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$OUTPUT_DIRECTORY" ] || [ -z "$REPLICATE_INDEX" ]; then
    echo "Error: INPUT_DIRECTORY, OUTPUT_DIRECTORY and REPLICATE_INDEX must be provided."
    show_help
    exit 1
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

Rscript /mnt/shiwei/bin/Replicate.R ${INPUT_DIRECTORY} ${OUTPUT_DIRECTORY} ${REPLICATE_INDEX}

echo "132 `get_time` Processing completed ..."


