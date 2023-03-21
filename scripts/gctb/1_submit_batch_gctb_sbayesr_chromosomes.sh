#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes:

# Before running SBayesR, filter SNPs in GWAS summary statistics.
# Although SBayesR will automatically filter SNPs to HapMap3 in LD matrix.


################################################################################
################################################################################
################################################################################



################################################################################
# Organize arguments.

path_file_gwas=${1} # full path and name to file for source GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_directory_ld_matrix=${2} # full path to parent directory for chromosome-specific Linkage Disequilibrium (LD) reference matrices in GCTB format
name_file_ld_matrix_prefix=${3} # prefix of name of file for chromosome-specific LD reference matrices
name_file_ld_matrix_suffix=${4} # suffix of name of file for chromosome-specific LD reference matrices
path_directory_product=${5} # full path to parent directory for product files from GCTB SBayesR
name_file_product_prefix=${6} # prefix of name of file for product files from GCTB SBayesR
name_file_product_suffix=${7} # suffix of name of file for product files from GCTB SBayesR
observations_variant=${8} # logical binary indicator of whether counts of observations are reliable and specific to each variant (SNP)
chromosome_x=${9} # whether to include Chromosome X
path_script_batch_run_sbayesr=${10} # full path to directory and file of script for execution of batch job for run of GCTB SBayesR
path_script_run_sbayesr=${11} # full path to directory and file of script for direct run of GCTB SBayesR
path_gctb=${12} # full path to directory and file for local executable installation of GCTB SBayesR
threads=${13} # count of concurrent or parallel process threads on node cores
report=${14} # whether to print reports



################################################################################
# Organize paths.

# Files.
path_file_batch_instances="${path_directory_product}/batch_instances.txt"
path_file_batch_out="${path_directory_product}/batch_out.txt"
path_file_batch_error="${path_directory_product}/batch_error.txt"

# Initialize directories and files.
mkdir -p $path_directory_product
rm $path_file_batch_instances
rm $path_file_batch_out
rm $path_file_batch_error



################################################################################
# Organize batch job instances.

# Iterate on relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define paths and names to files for chromosome-specific LD matrix.
  # Use base file name.
  path_file_base_ld_matrix="${path_directory_ld_matrix}/${name_file_ld_matrix_prefix}${chromosome}${name_file_ld_matrix_suffix}"
  # Define paths and names to files for results from SBayesR procedure.
  # Use base file name.
  path_file_base_product="${path_directory_product}/${name_file_product_prefix}${chromosome}${name_file_product_suffix}"

  # Define and append a new batch instance.
  instance="${chromosome};${path_file_base_ld_matrix};${path_file_base_product}"
  echo $instance >> $path_file_batch_instances
done

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "submit_batch_gctb_sbayesr_chromosomes.sh"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}
  echo "----------"
fi



################################################################################
# Submit batch to cluster scheduler for processing.
# Submit array of jobs in batch to Sun Grid Engine.
# Array batch indices must start at one (not zero).

if true; then
  qsub -t 1-${batch_instances_count}:1 \
  -o $path_file_batch_out \
  -e $path_file_batch_error \
  $path_script_batch_run_sbayesr \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_file_gwas \
  $observations_variant \
  $path_script_run_sbayesr \
  $path_gctb \
  $threads \
  $report
fi



#
