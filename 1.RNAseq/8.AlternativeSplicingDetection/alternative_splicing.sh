#! /usr/bin/bash



# Description : Detect the differentail splicing events by rMATs based on the bam files. 
# Author      : Zheting Zhang
# Time        : 2024-11-17



# Function to display help

display_help() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  --b1 <b1 file path>        Path to the b1 file containing BAM file paths seprated by comma"
  echo "  --b2 <b2 file path>        Path to the b2 file containing BAM file paths seprated by comma"
  echo "  --gtf <GTF file path>      Path to the GTF file"
  echo "  --od <output directory>    Output directory for the results"
  echo "  --help                     Display this help message"
  exit 0

}

# Parse command line arguments

while [[ $# -gt 0 ]]; do

  key="$1"
  case $key in

    --b1)
      b1_file="$2"
      shift # past argument

      shift # past value

      ;;
    --b2)
      b2_file="$2"
      shift # past argument

      shift # past value

      ;;
    --gtf)
      gtf_file="$2"
      shift # past argument

      shift # past value

      ;;
    --od)
      output_dir="$2"
      shift # past argument

      shift # past value

      ;;
    --help)
      display_help

      ;;
    *)
      echo "Error: Unknown option $1"
      display_help

      ;;
  esac

done

# Check if required arguments are provided

if [ -z "$b1_file" ] || [ -z "$b2_file" ] || [ -z "$gtf_file" ] || [ -z "$output_dir" ]; then

  echo "Error: All required arguments must be provided"
  display_help

  exit 1

fi

# Check if input files exist

if [ ! -f "$b1_file" ]; then

  echo "Error: b1 file $b1_file does not exist"
  exit 1

fi

if [ ! -f "$b2_file" ]; then

  echo "Error: b2 file $b2_file does not exist"
  exit 1

fi

if [ ! -f "$gtf_file" ]; then

  echo "Error: GTF file $gtf_file does not exist"
  exit 1

fi

# Run rmats.py command

mkdir "$output_dir"

mkdir "$output_dir"/tmp

rmats.py \
  --b1 "$b1_file" \
  --b2 "$b2_file" \
  --gtf "$gtf_file" \
  -t paired \
  --readLength 150 \
  --nthread 20 \
  --od "$output_dir" \
  --tmp "$output_dir"/tmp
  1>rmats.log 2>rmats.err &

echo "rmats.py command has started, output will be saved in $output_dir, log files are rmats.log and rmats.err"
