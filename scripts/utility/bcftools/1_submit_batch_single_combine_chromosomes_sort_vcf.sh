#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 06 June 2022
################################################################################
# Notes:
# This script collects source genotype files in Variant Call Format (VCF) format
# for chromosomes within a list and submits these for combination and sorting.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_directory_genotype_vcf_source=${1} # full path to directory with source genotype files in VCF format
prefix_file_genotype_vcf_source=${2} # file name prefix for source genotype file in VCF format
suffix_file_genotype_vcf_source=${3} # file name suffix for source genotype file in VCF format
chromosome_x=${4} # whether to include Chromosome X
path_file_list_source_vcf_files=${5} # full path to file with line-delimiter list of full paths to genotype files in VCF formats with BGZip compression and Tabix indices
path_file_genotype_vcf_product=${6} # full path to product genotype file in VCF format
threads=${7} # count of processing threads to use
path_promiscuity_scripts=${8} # full path to directory of general scripts
path_bcftools=${9} # full path to installation executable file of BCFTools
report=${10} # whether to print reports

###########################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_genotype_vcf_product .vcf.gz)"
path_directory_genotype_vcf_product="$(dirname $path_file_genotype_vcf_product)"
path_batch_instances="${path_directory_genotype_vcf_product}/batch_instances.txt"

# Scripts.
path_script_run_combine_sort="${path_promiscuity_scripts}/utility/bcftools/2_run_batch_single_combine_sort_vcf.sh"
path_script_combine_sort="${path_promiscuity_scripts}/utility/bcftools/3_combine_sort_vcf.sh"

# Initialize directory.
rm -r $path_directory_genotype_vcf_product
mkdir -p $path_directory_genotype_vcf_product
rm $path_file_list_source_vcf_files
rm $path_batch_instances

###########################################################################
# Collect paths to files for genotypes in VCF format.
# Iterate on source genotype files in VCF format for chromosomes.

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("21" "22" "X")
  #chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Define file names for chromosome.
  name_file_vcf_source_chromosome="${prefix_file_genotype_vcf_source}${chromosome}${suffix_file_genotype_vcf_source}"
  # Define full file paths for chromosome.
  path_file_vcf_source_chromosome="${path_directory_genotype_vcf_source}/${name_file_vcf_source_chromosome}"
  # Append file path to list.
  echo $path_file_vcf_source_chromosome >> $path_file_list_source_vcf_files
done

###########################################################################
# Define parameters for batch instances.

# Define and append a new batch instance.
instance="${path_file_list_source_vcf_files};${path_file_genotype_vcf_product}"
echo $instance >> $path_batch_instances

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_submit_batch_single_combine_chromosomes_sort_vcf.sh"
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
  "${path_script_run_combine_sort}" \
  $path_batch_instances \
  $batch_instances_count \
  $threads \
  $path_script_combine_sort \
  $path_bcftools \
  $report
fi



#
