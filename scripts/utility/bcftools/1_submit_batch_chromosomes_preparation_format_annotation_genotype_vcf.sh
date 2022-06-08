#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 26 May 2022
################################################################################
# Notes:
# This script finds within a parent directory all child genotype files in
# Variant Call Format (VCF).
# For each of these child genotype files in VCF format, the script calls
# another script that adjusts format and introduces annotations in BCFTools.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_directory_genotype_vcf_source=${1} # full path to directory with source genotype files in VCF format
prefix_file_genotype_vcf_source=${2} # file name prefix for source genotype file in VCF format
suffix_file_genotype_vcf_source=${3} # file name suffix for source genotype file in VCF format
chromosome_x=${4} # whether to include Chromosome X
path_directory_genotype_vcf_product=${5} # full path to directory for product genotype files in BCF format
prefix_file_genotype_vcf_product=${6} # file name prefix for product genotype file in VCF format
suffix_file_genotype_vcf_product=${7} # file name suffix for product genotype file in VCF format
path_file_genome_assembly_sequence=${8} # full path to genome assembly sequence file in FASTA format without compression that matches source and product genome assembly
path_file_chromosome_translations=${9} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_file_dbsnp_reference=${10} # full path to file for dbSNP reference in VCF format
threads=${11} # count of processing threads to use
path_promiscuity_scripts=${12} # full path to directory of general scripts
path_bcftools=${13} # full path to installation executable file of BCFTools
report=${14} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_directory_genotype_vcf_product}/batch_instances.txt"

# Scripts.
path_script_run_preparation_format_annotation="${path_promiscuity_scripts}/utility/bcftools/2_run_batch_preparation_format_annotation_vcf.sh"
path_script_drive_preparation_format_annotation="${path_promiscuity_scripts}/utility/bcftools/3_drive_prepare_format_annotate_vcf.sh"
#path_script_decompose_align_unique_sort="${path_promiscuity_scripts}/utility/bcftools/3_decompose_align_unique_sort_vcf.sh"
path_script_decompose_align_unique_sort="${path_promiscuity_scripts}/utility/bcftools/3_align_unique_sort_vcf.sh"
path_script_translate_chromosomes="${path_promiscuity_scripts}/utility/bcftools/3_translate_chromosome_identifiers_vcf.sh"
path_script_introduce_dbsnp_rsid="${path_promiscuity_scripts}/utility/bcftools/3_introduce_dbsnp_rsid_vcf.sh"

# Initialize directory.
rm -r $path_directory_genotype_vcf_product
mkdir -p $path_directory_genotype_vcf_product
rm $path_batch_instances

###########################################################################
# Iterate on source genotype files in VCF format for chromosomes.

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  #chromosomes=("1")
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Chromosome identifier for product file name.
  chromosome_lower="${chromosome,,}"
  # Define file names for chromosome.
  name_file_vcf_source_chromosome="${prefix_file_genotype_vcf_source}${chromosome}${suffix_file_genotype_vcf_source}"
  name_file_vcf_product_chromosome="${prefix_file_genotype_vcf_product}${chromosome_lower}${suffix_file_genotype_vcf_product}"
  # Define full file paths for chromosome.
  path_file_vcf_source_chromosome="${path_directory_genotype_vcf_source}/${name_file_vcf_source_chromosome}"
  path_file_vcf_product_chromosome="${path_directory_genotype_vcf_product}/${name_file_vcf_product_chromosome}"
  # Define and append a new batch instance.
  instance="${path_file_vcf_source_chromosome};${path_file_vcf_product_chromosome}"
  echo $instance >> $path_batch_instances
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_submit_batch_chromosomes_preparation_format_annotation_genotype_vcf.sh"
  echo "----------"
fi

################################################################################
# Submit batch instances to cluster scheduler.

# Read batch instances.
readarray -t batch_instances < $path_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}

# Execute batch with grid scheduler.
if true; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_directory_genotype_vcf_product}/batch_out.txt" -e "${path_directory_genotype_vcf_product}/batch_error.txt" \
  "${path_script_run_preparation_format_annotation}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_file_genome_assembly_sequence \
  $path_file_chromosome_translations \
  $path_file_dbsnp_reference \
  $threads \
  $path_script_drive_preparation_format_annotation \
  $path_script_decompose_align_unique_sort \
  $path_script_translate_chromosomes \
  $path_script_introduce_dbsnp_rsid \
  $path_bcftools \
  $report
fi



#
