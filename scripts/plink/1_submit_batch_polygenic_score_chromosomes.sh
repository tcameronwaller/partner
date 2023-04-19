#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 March 2023
# Date, last execution: 5 April 2023
# Review: TCW; 5 April 2023
################################################################################
# Note

# When calculating polygenic scores on large cohorts of genotypes, this process
# requires considerable computational resources and multiple hours of time.

# Write a batch submission script to handle the parallelization across
# chromosomes.

# Each set of SBayesR effects corresponds to one batch of jobs for each
# autosomal chromosome.



################################################################################
# Organize arguments.

path_file_source_effects=${1} # full path to source file in standard format of allelic effects across SNPs
path_directory_source_genotypes=${2} # full path to source directory for target genotypes in Variant Call Format (VCF)
name_file_genotypes_prefix=${3} # prefix in names of files for target genotypes in Variant Call Format (VCF)
name_file_genotypes_suffix=${4} # suffix in names of files for target genotypes in Variant Call Format (VCF)
path_directory_product=${5} # full path to parent directory for product files
name_file_product_prefix=${6} # prefix in name of product file for polygenic scores of allelic effects across target genotypes
name_file_product_suffix=${7} # suffix in name of product file for polygenic scores of allelic effects across target genotypes
report=${8} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_logs="${path_directory_product}/logs"

# Files.
path_file_batch_instances="${path_directory_product}/batch_instances.txt"
#path_file_batch_out="${path_directory_product}/batch_out.txt"
#path_file_batch_error="${path_directory_product}/batch_error.txt"

# Scripts.
path_script_run_batch="${path_directory_process}/promiscuity/scripts/plink/2_run_batch_polygenic_score.sh"
path_script_calculate_score="${path_directory_process}/promiscuity/scripts/plink/calculate_polygenic_score.sh"

# Initialize directories.
mkdir -p $path_directory_product
mkdir -p $path_directory_logs
cd $path_directory_product

# Initialize files.
rm $path_file_batch_instances
#rm $path_file_batch_out
#rm $path_file_batch_error

###########################################################################
# Organize parameters.

chromosome_x="false"
threads=4

###########################################################################
# Execute procedure.

# Iterate on relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
  #chromosomes=("14" "15" "16" "17" "18" "19" "20" "21" "22") # temporarily finish chromosomes
fi
for chromosome in "${chromosomes[@]}"; do
  # Organize paths and parameters.
  path_file_source_genotypes="${path_directory_source_genotypes}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
  name_base_file_product="${name_file_product_prefix}${chromosome}${name_file_product_suffix}"
  # Define parameters in array instance for batch job.
  instance="$path_file_source_effects;$path_file_source_genotypes;$name_base_file_product"
  echo $instance >> $path_file_batch_instances
done

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
index_array_maximum=$(batch_instances_count - 1)



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "1_submit_batch_polygenic_scores_chromosomes.sh"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$index_array_maximum]}
  echo "----------"
fi



################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Oracle Sun Grid Engine.
# Indices in array of batch jobs start at zero.

if true; then
  sbatch --array 0-${index_array_maximum}:1 --chdir $path_directory_product \
  $path_script_run_batch \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_directory_product \
  $threads \
  $report \
  $path_script_calculate_score
fi



#
