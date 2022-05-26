#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Note

# BCFTools specializes in tasks on genotype files in the Variant Call Format
# (VCF) that describe Single Nucleotide Polymorphisms (SNPs) across samples.

# BCFTools documentation
# https://samtools.github.io/bcftools/bcftools.html

# The "concat" command in BCFTools requires that the multiple VCF files all have
# the same samples.

################################################################################

# TODO: TCW; 26 May 2022
# TODO: the preparation of genotype files for combination could be split into an array batch job
#


################################################################################
# Organize arguments.

path_directory_source_genotype_vcf=${1} # full path to directory with source genotype files in VCF format
prefix_file_source_genotype_vcf=${2} # file name prefix for source genotype file in VCF format
suffix_file_source_genotype_vcf=${3} # file name suffix for source genotype file in VCF format
chromosome_x=${4} # whether to include Chromosome X
path_directory_product_genotype_vcf=${5}
threads=${6} # count of processing threads to use
path_promiscuity_scripts=${7} # full path to directory of general scripts
path_bcftools=${8} # full path to installation executable file of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_directory_product_genotype_vcf_temporary="${path_directory_product_genotype_vcf}/temporary"
path_file_list_files_combination="${path_directory_product_genotype_vcf_temporary}/list_files_combination.txt"

# Scripts.
path_script_preparation="${path_promiscuity_scripts}/utility/bcftools/2_convert_vcf_bcf_remove_duplicates_sort_samples_records.sh"
path_script_combination="${path_promiscuity_scripts}/utility/bcftools/3_combine_sort_multiple_bcf_convert_vcf.sh"

###########################################################################
# Iterate on source genotype files in VCF format for chromosomes.

# Initialize directory.
rm -r $path_directory_product_genotype_vcf
mkdir -p $path_directory_product_genotype_vcf
rm -r $path_directory_product_genotype_vcf_temporary
mkdir -p $path_directory_product_genotype_vcf_temporary

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
# Iterate on chromosomes.
for chromosome in "${chromosomes[@]}"; do
  # Define file names for chromosome.
  name_file_source_vcf_chromosome="${prefix_file_source_genotype_vcf}${chromosome}${suffix_file_source_genotype_vcf}"
  name_file_intermediate_bcf_chromosome="genotype_chromosome${chromosome}_bcf"
  name_file_intermediate_remove_duplicates_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates"
  name_file_intermediate_list_samples_chromosome="list_samples_chromosome${chromosome}.txt"
  name_file_intermediate_sort_samples_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates_sort_samples"
  name_file_intermediate_sort_records_chromosome="genotype_chromosome${chromosome}_bcf_remove_duplicates_sort_samples_records"
  # Define full file paths for chromosome.
  path_file_source_vcf_chromosome="${path_directory_source_genotype_vcf}/${name_file_source_vcf_chromosome}"
  path_file_intermediate_bcf_chromosome="${path_directory_product_genotype_vcf_temporary}/${name_file_intermediate_bcf_chromosome}"
  path_file_intermediate_remove_duplicates_chromosome="${path_directory_product_genotype_vcf_temporary}/${name_file_intermediate_remove_duplicates_chromosome}"
  path_file_intermediate_list_samples_chromosome="${path_directory_product_genotype_vcf_temporary}/${name_file_intermediate_list_samples_chromosome}"
  path_file_intermediate_sort_samples_chromosome="${path_directory_product_genotype_vcf_temporary}/${name_file_intermediate_sort_samples_chromosome}"
  path_file_intermediate_sort_records_chromosome="${path_directory_product_genotype_vcf_temporary}/${name_file_intermediate_sort_records_chromosome}"
  # Call script for preparation of the genotype file.
  /usr/bin/bash "${path_script_preparation}" \
  $path_file_source_vcf_chromosome \
  $path_directory_product_genotype_vcf_temporary \
  $path_file_intermediate_bcf_chromosome \
  $path_file_intermediate_remove_duplicates_chromosome \
  $path_file_intermediate_list_samples_chromosome \
  $path_file_intermediate_sort_samples_chromosome \
  $path_file_intermediate_sort_records_chromosome \
  $threads \
  $path_bcftools \
  $report
  # Define and append a new batch instance.
  instance="${path_file_intermediate_sort_records_chromosome}"
  echo $instance >> $path_file_list_files_combination
done


# 5. combine BCF files
# 6. sort SNP records
# 8. convert back to VCF format


# Call script for combination of genotype files.
#/usr/bin/bash "${path_script_combination}" \
#$path_file_list_files_combination \
#$threads \
#$path_bcftools \
#$report

# Remove temporary, intermediate files.
#rm -r $path_directory_product_genotype_vcf_temporary



#
