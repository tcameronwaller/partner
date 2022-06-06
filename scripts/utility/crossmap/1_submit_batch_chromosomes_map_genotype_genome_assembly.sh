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

# TODO: TCW; 26 May 2022
# TODO: 1. find genotype VCF files for all chromosomes
# TODO: 2. combine these together into a single VCF file
# TODO: 3. sort that single VCF file
# TODO: 4. perform assembly translation in CrossMap
# TODO: 5. sort the VCF file after translation
# TODO: 6. split by chromosome


################################################################################
# Organize arguments.

path_file_source_vcf=${1} # full path to source file in VCF format
path_file_product_vcf=${2} # full path to product file in VCF format
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_genome_assembly_sequence=${4} # full path to genome assembly sequence file in FASTA format without compression that matches product genome assembly
threads=${5} # count of processing threads to use
path_promiscuity_scripts=${6} # full path to directory of general scripts
path_environment_crossmap=${7} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${8} # full path to installation executable file of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_directory_product="$(dirname $path_file_product_vcf)"
path_batch_instances="${path_directory_product}/batch_instances.txt"

# Scripts.
path_script_run_map_assembly="${path_promiscuity_scripts}/utility/crossmap/2_run_batch_map_genotype_genome_assembly.sh"
path_script_map_assembly="${path_promiscuity_scripts}/utility/crossmap/3_map_genotype_genome_assembly_vcf.sh"

# Initialize directory.
rm -r $path_directory_product
mkdir -p $path_directory_product
rm -r $path_batch_instances

###########################################################################
# Define parameters for batch instances.

# Define and append a new batch instance.
instance="${path_file_source_vcf};${path_file_product_vcf}"
echo $instance >> $path_batch_instances

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_submit_batch_single_translate_genome_assembly.sh"
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
  "${path_script_run_map_assembly}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_assembly_translation_chain \
  $path_genome_assembly_sequence \
  $threads \
  $path_script_map_assembly \
  $path_environment_crossmap \
  $path_bcftools \
  $report
fi



#
