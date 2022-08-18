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

path_file_gwas_summary=${1} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${2} # integer count of samples for the GWAS study
path_directory_genotypes_snp_bim=${3} # full path to directory for relevant SNPs in BIM format
prefix_file_genotype_snp_bim=${4} # file name prefix for list of relevant target SNPs in BIM format without '.bim' suffix
suffix_file_genotype_snp_bim=${5} # file name suffix for list of relevant target SNPs in BIM format without '.bim' suffix
path_directory_genetic_reference_prscs=${6} # full path to directory for genetic references
population_ancestry=${7} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
path_directory_allele_effect=${8} # full path to directory for product reports on posterior allele effect size estimates
name_file_allele_effect_prefix=${9} # prefix for name of product report file
chromosome_x=${10} # whether to include Chromosome X
threads=${11} # count of processing threads to use
path_directory_promiscuity_scripts=${12} # full path to script for estimation of allelic posterior effects in PRS-CSX
path_environment_prscs=${13} # full path to Python 3 environment with installation of CrossMap
path_prscsx=${14} # full path to installation executable file of PRS-CSX
report=${15} # whether to print reports

###########################################################################
# Organize paths.

path_file_batch_instances="${path_directory_allele_effect}/batch_instances.txt"

# Scripts.
path_script_run_prscs_estimate="${path_directory_promiscuity_scripts}/utility/prscs_polygenic_score/2_run_batch_prscs_estimate_allelic_effects.sh"
path_script_prscs_estimate_allelic_effects="${path_directory_promiscuity_scripts}/utility/prscs_polygenic_score/3_estimate_prscs_allelic_effects.sh"

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
  name_bim_file_chromosome="${prefix_file_genotype_snp_bim}${chromosome}${suffix_file_genotype_snp_bim}"
  path_file_genotype_snp_bim_prefix="${path_directory_genotypes_snp_bim}/${name_bim_file_chromosome}"
  # Define name of product file for chromosome.
  name_file_allele_effect="${name_file_allele_effect_prefix}_chromosome${chromosome}"
  # Define and append a new batch instance.
  instance="${chromosome};${path_file_genotype_snp_bim_prefix};${name_file_allele_effect}"
  echo $instance >> $path_file_batch_instances
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
readarray -t batch_instances < $path_file_batch_instances
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
  "${path_directory_allele_effect}/batch_out.txt" -e "${path_directory_allele_effect}/batch_error.txt" \
  "${path_script_run_prscs_estimate}" \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_file_gwas_summary \
  $count_gwas_samples \
  $path_directory_genetic_reference_prscs \
  $population_ancestry \
  $path_directory_allele_effect \
  $threads \
  $path_script_prscs_estimate_allelic_effects \
  $path_environment_prscs \
  $path_prscsx \
  $report
fi



#
