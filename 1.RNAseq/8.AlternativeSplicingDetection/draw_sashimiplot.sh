#! /usr/bin/bash

# Description : Draw the sashimi plot to visualize the differential splicing events based on the bam files.
# Author      : Zheting Zhang
# Time        : 2024-11-17




# Function to display help

display_help() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  --b1 <b1 file path>        Path to the b1 file containing BAM file paths seprated by comma"
  echo "  --b2 <b2 file path>        Path to the b2 file containing BAM file paths seprated by comma"
  echo "  --gff <GFF file path>      Path to the GFF file"
  echo "  --od <output directory>    Output directory for the results"
  echo "  --r <interval>             Interval of interests,chromosome:orientation:start:end. example: chr11:-:4752513:4854957"
  echo "  --g <group file path>           Path to the file containing group information of samples. example:  WT: 1-3
                                                                                                      Mut: 4-6    "	
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
    --gff)
      gff_file="$2"
      shift # past argument

      shift # past value

      ;;
    --od)
      output_dir="$2"
      shift # past argument

      shift # past value

      ;;
    --r)
      interval="$2"
      shift # past argument

      shift # past value

      ;;
    --g) 
      groupinfo="$2"
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

if [ -z "$b1_file" ] || [ -z "$b2_file" ] || [ -z "$gff_file" ] || [ -z "$output_dir" ]||[ -z "$interval" ]||[ -z "$groupinfo"]; then

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

if [ ! -f "$gff_file" ]; then

  echo "Error: GFF file $gff_file does not exist"
  exit 1

fi

if [ ! -f "$groupinfo" ]; then

  echo "Error: group infomation file $groupinfo does not exist"
  exit 1

fi

# Run rmats2sashimiplot command

mkdir "$output_dir";

mkdir "$output_dir"/tmp;

b1info=$(cat "$b1_file");

b2info=$(cat "$b2_file");

sample1=$(cut -f1 -d " " "$groupinfo"|sed -e "s/://g"|head -n1);

sample2=$(cut -f1 -d " " "$groupinfo"|sed -e "s/://g"|tail -n1);

rmats2sashimiplot \
  --b1 "$b1info" \
  --b2 "$b2info" \
  -c "$interval":"$gff_file" \
  --intron_s 5 \
  -o "$output_dir" \
  --group-info "$groupinfo" \
  --l1 "$sample1" --l2 "$sample2"

