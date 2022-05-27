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

################################################################################
# Organize arguments.
path_file_source_vcf_chromosome=${1} # full path to source genotype file in VCF format
path_directory_product_temporary_chromosome=${2} # full path to directory for temporary intermediate files
path_file_intermediate_format_chromosome=${3} # full path to intermediate file
path_file_intermediate_remove_duplicates_chromosome=${4} # full path to intermediate file
path_file_intermediate_list_samples_chromosome=${5} # full path to intermediate file
path_file_intermediate_sort_samples_chromosome=${6} # full path to intermediate file
path_file_intermediate_sort_records_chromosome=${7} # full path to intermediate file
path_file_product_bcf_chromosome=${8} # full path to product genotype file in BCF format
threads=${9} # count of processing threads to use
path_bcftools=${10} # full path to installation executable file of BCFTools
report=${11} # whether to print reports

################################################################################
# Organize paths.

# Initialize directory.
rm -r $path_directory_product_temporary_chromosome
mkdir -p $path_directory_product_temporary_chromosome

################################################################################
# Prepare genotype files for combination.



# Convert genotype files from VCF format with BGZip compression to working
# format.
# The BCF format allows for greater performance in BCFTools.
# Read file in VCF format with BGZip compression in BCFTools.
# Write file in BCF format without compression or other format.
echo "----------"
echo "$path_file_source_vcf_chromosome"
echo "----------"
$path_bcftools \
view \
--output $path_file_intermediate_format_chromosome \
--output-type v \
--threads $threads \
$path_file_source_vcf_chromosome
echo "----------"
echo "$path_file_intermediate_format_chromosome"
echo "----------"


# Remove duplicate records for SNPs or other genetic features.
# I think that option "--rm-dup exact" evokes similar logic to option
# "--collapse none".
# I think both of these options require exact match of chromosome, position, and
# all alleles in order to consider records redundant.
#$path_bcftools \
#norm \
#--rm-dup exact \
#--output $path_file_intermediate_remove_duplicates_chromosome \
#--output-type u \
#--threads $threads \
#$path_file_intermediate_bcf_chromosome

# Sort samples.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf
# Extract sample identifiers from genotype file and sort these externally.
$path_bcftools \
query \
--list-samples \
$path_file_intermediate_format_chromosome | sort > $path_file_intermediate_list_samples_chromosome
# Sort samples within genotype file.
$path_bcftools \
view \
--samples-file $path_file_intermediate_list_samples_chromosome \
--output $path_file_intermediate_sort_samples_chromosome \
--output-type v \
--threads $threads \
$path_file_intermediate_format_chromosome
echo "----------"
echo "$path_file_intermediate_sort_samples_chromosome"
echo "----------"
# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 5G \
--output $path_file_intermediate_sort_records_chromosome \
--output-type v \
--temp-dir $path_directory_product_genotype_vcf_temporary \
$path_file_intermediate_sort_samples_chromosome
echo "----------"
echo "$path_file_intermediate_sort_records_chromosome"
echo "----------"

# Copy to more permanent product file.
cp $path_file_intermediate_sort_records_chromosome $path_file_product_bcf_chromosome

# Remove temporary, intermediate files.
#rm -r $path_directory_product_temporary_chromosome



#
