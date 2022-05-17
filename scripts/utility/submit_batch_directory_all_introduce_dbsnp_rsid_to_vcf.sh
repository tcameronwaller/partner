#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Notes:
# This script finds within a parent directory all child genotype files in
# Variant Call Format (VCF).
# For each of these child genotype files in VCF format, the script calls
# another script to extract information about Single Nucleotide polymorphisms
# (SNPs) to a new file in PLINK2 BIM format.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_genotype_source_vcf_container=${1} # full path to parent directory with source genotype files in VCF format
pattern_genotype_source_vcf_file=${2} # string glob pattern by which to recognize source genotype files in VCF format
path_genotype_product_vcf_container=${3} # full path to parent directory for product genotype files in BIM format
path_promiscuity_scripts=${4} # full path to directory of general scripts
report=${5} # whether to print reports

###########################################################################
# Organize paths.

path_batch_instances="${path_genotype_product_vcf_container}/batch_instances.txt"

# Scripts.
path_script_run_dbsnp_rsid_to_vcf="${path_promiscuity_scripts}/utility/run_batch_introduce_dbsnp_rsid_to_vcf.sh"
path_script_dbsnp_rsid_to_vcf="${path_promiscuity_scripts}/utility/introduce_dbsnp_rsid_to_vcf.sh"

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
    # TODO: append a new batch instance...
    instance="${path_vcf_source};${path_vcf_product}"
    echo $instance >> $path_batch_instances

    # Convert information from genotype files in VCF format to BIM format.
    /usr/bin/bash "${path_script_dbsnp_rsid_to_vcf}" \
    $path_genotype_source_vcf \
    $path_genotype_product_bim_container \
    $name_file_product_bim \
    $report
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
  "${path_script_run_dbsnp_rsid_to_vcf}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_script_dbsnp_rsid_to_vcf
fi
