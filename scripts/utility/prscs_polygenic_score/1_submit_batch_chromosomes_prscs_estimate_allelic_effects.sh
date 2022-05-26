#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 24 May 2022
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


# TODO: define "name_file_product" with chromosome number
prefix_genotype_snp_bim_file=${5} # full path to file name prefix for list of relevant target SNPs in BIM format without '.bim' suffix
suffix_genotype_snp_bim_file=${6} # full path to file name prefix for list of relevant target SNPs in BIM format without '.bim' suffix
name_file_product_prefix=${10} # prefix for name of product report file
chromosome_x=${9} # whether to collect GWAS summary statistics report for Chromosome X
path_genotype_snp_bim_directory=${0} # full path to
path_promiscuity_scripts=${7} # full path to directory of general scripts

###########################################################################
# Organize paths.

path_batch_instances="${path_genotype_product_vcf_container}/batch_instances.txt"

# Scripts.
path_script_run_prscs_estimate="${path_promiscuity_scripts}/utility/prscs_polygenic_score/2_run_batch_prscs_estimate_allelic_effects.sh"
path_script_prscs_estimate_allelic_effects="${path_promiscuity_scripts}/utility/prscs_polygenic_score/estimate_prscs_allelic_effects.sh"

###########################################################################
# Find source genotype files in VCF format within container directory.

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define paths to SNP BIM reference for chromosome.
  name_bim_file_chromosome="${prefix_genotype_snp_bim_file}${chromosome}${suffix_genotype_snp_bim_file}"
  path_snp_relevance_bim_prefix="${path_genotype_snp_bim_directory}/${name_bim_file_chromosome}"
  # Define name of product file for chromosome.
  name_file_product="${name_file_product_prefix}_chromosome${chromosome}"
  # Define and append a new batch instance.
  instance="${chromosome};${path_snp_relevance_bim_prefix};${name_file_product}"
  echo $instance >> $path_batch_instances
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
  "${path_genotype_product_vcf_container}/batch_out.txt" -e "${path_genotype_product_vcf_container}/batch_error.txt" \
  "${path_script_run_prscs_estimate}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_source_gwas_summary \
  $count_gwas_samples \
  $path_genetic_reference_prscs \
  $population_ancestry \
  $path_product_allele_effect_directory \
  $threads \
  $path_script_prscs_estimate_allelic_effects \
  $path_environment_prscs \
  $path_prscsx \
  $report
fi



#
