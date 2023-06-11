#!/bin/bash
# Author: WANG Xuedong
# Date: 2023-06-10
# Description: A SNP data prepare script for the MHC project

# Global variables
readonly SCRIPT_NAME="$(basename "$0")"
readonly LOG_FILE="./${SCRIPT_NAME%.sh}.log"

# Function definitions

# Usage: display_usage
# Description: Display script usage information
display_usage() {
  echo "Usage: $SCRIPT_NAME [options] <save_path>"
  echo
  echo "Options:"
  echo "  -h, --help    Display this help message and exit"
  echo
  echo "Arguments:"
  echo "  save_path          save path"
}

# Usage: log MESSAGE
# Description: Log the given MESSAGE to a log file
log() {
  local message="$1"
  echo "$(date): $message" >> "$LOG_FILE"
}

# Main script logic

# Check if the user requested help
if [[ "$#" -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  display_usage
  exit 0
fi

# Check if the correct number of arguments is provided
if [[ "$#" -ne 1 ]]; then
  echo "Error: Invalid number of arguments"
  echo
  display_usage
  exit 1
fi

save_path="$1"

# Perform the main script logic
log "Starting script with arguments: $save_path"
echo "Processing arguments $save_path..."



# download 1000 genomes data
cd $save_path
for i in {1..22}
do
    echo "Downloading chr$i..."
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz ./
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz.tbi ./
    echo "finish downloading chr$i"
done

# extract MHC region
chr6_vcf=$save_path/CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.filtered.shapeit2-duohmm-phased.vcf.gz
MHC_vcf=$save_path/MHC.vcf.gz
bcftools view -r chr6:29720000-33130000 $chr6_vcf -Oz -o $MHC_vcf
bcftools index $MHC_vcf
# split other chromosomes by MHC region snp number
snp_cnt=$(bcftools view -H $MHC_vcf | wc -l)
echo "MHC region has $snp_cnt SNPs"
for i in {1..22}
do
    split_dir=split_chr$i
    mkdir -p $split_dir
    cd $split_dir
    bcftools query -f'%CHROM\t%POS\n' ../CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz | split -l $sp_cnt
    for file in x*; do bcftools view -T $file -Oz ../CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz > $file.vcf.gz; done
    cd ../
done

echo "finish jobs"