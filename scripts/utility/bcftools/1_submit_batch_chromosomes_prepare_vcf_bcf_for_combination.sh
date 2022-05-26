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

path_directory_source_genotype_vcf=${1} # full path to directory with source genotype files in VCF format
prefix_file_source_genotype_vcf=${2} # file name prefix for source genotype file in VCF format
suffix_file_source_genotype_vcf=${3} # file name suffix for source genotype file in VCF format
chromosome_x=${4} # whether to include Chromosome X
path_directory_product_genotype_bcf=${5} # full path to directory for product genotype files in BCF format
threads=${6} # count of processing threads to use
path_promiscuity_scripts=${7} # full path to directory of general scripts
path_bcftools=${8} # full path to installation executable file of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_directory_product_genotype_bcf}/batch_instances.txt"
path_file_list_files_combination="${path_directory_product_genotype_bcf}/list_files_chromosomes_combination.txt"

# Scripts.
path_script_run_preparation="${path_promiscuity_scripts}/utility/bcftools/2_run_batch_chromosome_prepare_vcf_bcf_for_combination.sh"
path_script_preparation="${path_promiscuity_scripts}/utility/bcftools/3_convert_vcf_bcf_remove_duplicates_sort_samples_records.sh"

###########################################################################
# Iterate on source genotype files in VCF format for chromosomes.

# Initialize directory.
rm -r $path_directory_product_genotype_vcf
mkdir -p $path_directory_product_genotype_vcf
rm -r $path_directory_product_genotype_vcf_temporary
mkdir -p $path_directory_product_genotype_vcf_temporary

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Define full path to temporary directory for the chromosome.
  path_directory_product_temporary_chromosome="${path_directory_product_genotype_bcf}/temporary_chromosome${chromosome}"
  # Define file names for chromosome.
  name_file_source_vcf_chromosome="${prefix_file_source_genotype_vcf}${chromosome}${suffix_file_source_genotype_vcf}"
  name_file_intermediate_bcf_chromosome="genotype_chromosome${chromosome}_bcf"
  name_file_intermediate_remove_duplicates_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates"
  name_file_intermediate_list_samples_chromosome="list_samples_chromosome${chromosome}.txt"
  name_file_intermediate_sort_samples_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates_sort_samples"
  name_file_intermediate_sort_records_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates_sort_samples_records"
  name_file_product_bcf_chromosome="genotype_chromosome${chromosome}_bcf.bcf"
  # Define full file paths for chromosome.
  path_file_source_vcf_chromosome="${path_directory_source_genotype_vcf}/${name_file_source_vcf_chromosome}"
  path_file_intermediate_bcf_chromosome="${path_directory_product_temporary_chromosome}/${name_file_intermediate_bcf_chromosome}"
  path_file_intermediate_remove_duplicates_chromosome="${path_directory_product_temporary_chromosome}/${name_file_intermediate_remove_duplicates_chromosome}"
  path_file_intermediate_list_samples_chromosome="${path_directory_product_temporary_chromosome}/${name_file_intermediate_list_samples_chromosome}"
  path_file_intermediate_sort_samples_chromosome="${path_directory_product_temporary_chromosome}/${name_file_intermediate_sort_samples_chromosome}"
  path_file_intermediate_sort_records_chromosome="${path_directory_product_temporary_chromosome}/${name_file_intermediate_sort_records_chromosome}"
  path_file_product_bcf_chromosome="${path_directory_product_genotype_bcf}/${name_file_product_bcf_chromosome}"
  # Define and append a new batch instance.
  instance="${chromosome};${path_file_source_vcf_chromosome};${path_directory_product_temporary_chromosome};${path_file_intermediate_bcf_chromosome};${path_file_intermediate_remove_duplicates_chromosome};${path_file_intermediate_list_samples_chromosome};${path_file_intermediate_sort_samples_chromosome};${path_file_intermediate_sort_records_chromosome};${path_file_product_bcf_chromosome}"
  echo $instance >> $path_batch_instances
  # Collect list of genotype files in BCF format for subsequent combination.
  instance="${path_file_product_bcf_chromosome}"
  echo $instance >> $path_file_list_files_combination
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "submit_batch_chromosomes_prscs_estimate_allelic_effects.sh"
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
  "${path_directory_product_genotype_bcf}/batch_out.txt" -e "${path_directory_product_genotype_bcf}/batch_error.txt" \
  "${path_script_run_preparation}" \
  $path_batch_instances \
  $batch_instances_count \
  $threads \
  $path_script_preparation \
  $path_bcftools \
  $report
fi



#
