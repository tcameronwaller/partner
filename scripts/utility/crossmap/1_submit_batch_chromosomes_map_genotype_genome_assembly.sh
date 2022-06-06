#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 06 June 2022
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
path_file_assembly_translation_chain=${8} # full path to chain file for assembly translation
path_file_genome_assembly_sequence=${9} # full path to genome assembly sequence file in FASTA format without compression that matches product genome assembly
threads=${10} # count of processing threads to use
path_promiscuity_scripts=${11} # full path to directory of general scripts
path_environment_crossmap=${12} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${13} # full path to installation executable file of BCFTools
report=${14} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_directory_genotype_vcf_product}/batch_instances.txt"

# Scripts.
path_script_run_map_assembly="${path_promiscuity_scripts}/utility/crossmap/2_run_batch_map_genotype_genome_assembly.sh"
path_script_map_assembly="${path_promiscuity_scripts}/utility/crossmap/3_map_genotype_genome_assembly_vcf.sh"

# Initialize directory.
rm -r $path_directory_genotype_vcf_product
mkdir -p $path_directory_genotype_vcf_product
rm $path_batch_instances

###########################################################################
# Define parameters for batch instances.
# Iterate on source genotype files in VCF format for chromosomes.

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  #chromosomes=("21" "22" "X")
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Define file names for chromosome.
  name_file_vcf_source_chromosome="${prefix_file_genotype_vcf_source}${chromosome}${suffix_file_genotype_vcf_source}"
  name_file_vcf_product_chromosome="${prefix_file_genotype_vcf_product}${chromosome}${suffix_file_genotype_vcf_product}"
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
  echo "1_submit_batch_chromosomes_map_genotype_genome_assembly.sh"
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
  "${path_script_run_map_assembly}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_file_assembly_translation_chain \
  $path_file_genome_assembly_sequence \
  $threads \
  $path_script_map_assembly \
  $path_environment_crossmap \
  $path_bcftools \
  $report
fi



#
