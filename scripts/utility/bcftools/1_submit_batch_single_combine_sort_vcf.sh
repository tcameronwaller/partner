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

path_file_list_source_vcf_files=${1} # full path to file with line-delimiter list of full paths to genotype files in VCF formats with BGZip compression and Tabix indices
path_file_product_vcf=${2} # full path to product file in VCF format with BGZip compression
threads=${3} # count of processing threads to use
path_promiscuity_scripts=${4} # full path to directory of general scripts
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

###########################################################################
# Organize paths.

path_directory_product="$(dirname $path_file_product_vcf)"
path_batch_instances="${path_directory_product}/batch_instances.txt"

# Scripts.
path_script_run_combine_sort="${path_promiscuity_scripts}/utility/bcftools/2_run_batch_single_combine_sort_vcf.sh"
path_script_combine_sort="${path_promiscuity_scripts}/utility/bcftools/3_combine_sort_vcf.sh"

# Initialize directory.
#rm -r $path_directory_product
mkdir -p $path_directory_product
rm -r $path_batch_instances

###########################################################################
# Define parameters for batch instances.

# Define and append a new batch instance.
instance="${path_file_list_source_vcf_files};${path_file_product_vcf}"
echo $instance >> $path_batch_instances

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_submit_batch_single_combine_sort_vcf.sh"
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
  "${path_directory_product}/batch_out.txt" -e "${path_directory_product}/batch_error.txt" \
  "${path_script_run_combine_sort}" \
  $path_batch_instances \
  $batch_instances_count \
  $threads \
  $path_script_combine_sort \
  $path_bcftools \
  $report
fi



#
