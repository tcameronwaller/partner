#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 18 May 2022
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
path_genotype_source_vcf_container=${1} # full path to parent directory with source genotype files in VCF format
pattern_genotype_source_vcf_file=${2} # string glob pattern by which to recognize source genotype files in VCF format
path_genotype_product_vcf_container=${3} # full path to parent directory for product genotype files in BIM format
path_chromosome_translations=${4} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_dbsnp_reference=${5} # full path to file for dbSNP reference in VCF format
path_promiscuity_scripts=${6} # full path to directory of general scripts
threads=${7} # count of processing threads to use
path_bcftools=${8} # full path to installation of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_genotype_product_vcf_container}/batch_instances.txt"

# Scripts.
path_script_run_vcf_format_annotation="${path_promiscuity_scripts}/utility/bcftools/2_run_batch_vcf_format_annotation.sh"

###########################################################################
# Find source genotype files in VCF format within container directory.

cd $path_genotype_source_vcf_container
for path_file in `find . -maxdepth 1 -mindepth 1 -type f -name "$pattern_genotype_source_vcf_file"`; do
  if [[ -f "$path_file" ]]; then
    # Current content item is a file.
    # Extract directory's base name.
    name_file="$(basename -- $path_file)"
    echo $name_file
    path_vcf_source="${path_genotype_source_vcf_container}/${name_file}"
    path_vcf_product="${path_genotype_product_vcf_container}/${name_file}"
    # Define and append a new batch instance.
    instance="${path_vcf_source};${path_vcf_product}"
    echo $instance >> $path_batch_instances
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "submit_batch_directory_all_introduce_dbsnp_rsid_to_vcf.sh"
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
  "${path_script_run_vcf_format_annotation}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_chromosome_translations \
  $path_dbsnp_reference \
  $threads \
  $path_bcftools \
  $report
fi



#
