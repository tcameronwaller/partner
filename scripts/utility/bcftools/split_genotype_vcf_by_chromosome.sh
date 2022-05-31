#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 31 May 2022
################################################################################
# Notes:
# This script reads information in BCFTools from a single genotype file in
# Variant Call Format (VCF).
# This script extracts records for genetic features on each chromosome and
# writes these to separate files in VCF format.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_file_source_genotype_vcf=${1} # full path to directory with source genotype files in VCF format
chromosome_x=${2} # whether to include Chromosome X
path_directory_product_vcf=${3} # full path to directory for product genotype files in VCF format
prefix_file_product_vcf=${4} # file name suffix for product genotype files in VCF format
threads=${6} # count of processing threads to use
path_bcftools=${8} # full path to installation executable file of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_directory_product_genotype_bcf}/batch_instances.txt"
path_file_list_files_combination="${path_directory_product_genotype_bcf}/list_files_chromosomes_combination.txt"

# Scripts.
path_script_combine_sort_split="${path_promiscuity_scripts}/utility/bcftools/2_combine_bcf_sort_split_convert_vcf.sh"

###########################################################################
# Iterate on source genotype files in VCF format for chromosomes.

# Initialize directory.
rm -r $path_directory_product_genotype_bcf
mkdir -p $path_directory_product_genotype_bcf

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX")
else
  chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Define full path to the file for the chromosome.
  chromosome_simple="${chromosome//chr/""}"
  name_file_product_vcf_chromosome="${prefix_file_product_vcf}${chromosome_simple}.vcf.gz"
  path_file_product_vcf_chromosome="${path_directory_product_vcf}/${name_file_product_vcf_chromosome}"

  # Extract genotype records for genetic features on current chromosome.
  $path_bcftools \
  view \
  --regions $chromosome \
  --output $path_file_product_vcf_chromosome \
  --output-type z9 \
  --threads $threads \
  $path_file_source_genotype_vcf

  # Create Tabix index for file in VCF format with BGZip compression.
  # BCFTools is unable to create a Tabix index for files in BCF format.
  # Some commands in BCFTools require this Tabix index to read a file.
  $path_bcftools \
  index \
  --force \
  --tbi \
  --threads $threads \
  $path_file_product_vcf_chromosome
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "split_genotype_vcf_by_chromosome.sh"
  echo "----------"
fi



#
