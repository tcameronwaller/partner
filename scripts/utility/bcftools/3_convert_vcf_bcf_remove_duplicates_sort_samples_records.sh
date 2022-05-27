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
path_file_source_vcf=${1} # full path to source genotype file in VCF format
path_directory_product_temporary=${2} # full path to directory for temporary intermediate files, keeps chromosomes separate for parallel computing
path_file_product_bcf=${3} # full path to product genotype file in BCF format
threads=${4} # count of processing threads to use
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product_bcf .bcf.gz)"
path_directory_product="$(dirname $path_file_product_bcf)"
path_file_temporary_bcf="${path_directory_product_temporary}/${name_base_file_product}.bcf"
path_file_temporary_remove_replicates="${path_directory_product_temporary}/${name_base_file_product}_replicates.bcf"
path_file_temporary_list_samples="${path_directory_product_temporary}/${name_base_file_product}_list_samples.txt"
path_file_temporary_sort_samples="${path_directory_product_temporary}/${name_base_file_product}_sort_samples.bcf"
path_file_temporary_sort_records="${path_directory_product_temporary}/${name_base_file_product}_sort_records.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

################################################################################
# Prepare genotype files for combination.

# Convert genotype files from VCF format with BGZip compression to working
# format.
# The BCF format allows for greater performance in BCFTools.
# Read file in VCF format with BGZip compression in BCFTools.
# Write file in BCF format without compression or other format.
echo "----------"
echo "$path_file_source_vcf"
echo "----------"
$path_bcftools \
view \
--output $path_file_temporary_bcf \
--output-type u \
--threads $threads \
$path_file_source_vcf
echo "----------"
echo "$path_file_temporary_bcf"
echo "----------"

# Remove duplicate records for SNPs or other genetic features.
# I think that option "--rm-dup exact" evokes similar logic to option
# "--collapse none".
# I think both of these options require exact match of chromosome, position, and
# all alleles in order to consider records redundant.
$path_bcftools \
norm \
--rm-dup exact \
--output $path_file_temporary_remove_replicates \
--output-type u \
--threads $threads \
$path_file_temporary_bcf
echo "----------"
echo "$path_file_temporary_remove_replicates"
echo "----------"

# Sort samples.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf
# Extract sample identifiers from genotype file and sort these externally.
$path_bcftools \
query \
--list-samples \
$path_file_temporary_remove_replicates | sort > $path_file_temporary_list_samples
# Sort samples within genotype file.
$path_bcftools \
view \
--samples-file $path_file_temporary_list_samples \
--output $path_file_temporary_sort_samples \
--output-type u \
--threads $threads \
$path_file_temporary_remove_replicates
echo "----------"
echo "$path_file_temporary_sort_samples"
echo "----------"

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 4G \
--output $path_file_temporary_sort_records \
--output-type u \
--temp-dir $path_directory_product_temporary \
$path_file_temporary_sort_samples
echo "----------"
echo "$path_file_temporary_sort_records"
echo "----------"

# Convert genotype files from BCF format without compression to BCF format with
# BGZip compression.
# The BCF format allows for greater performance in BCFTools.
$path_bcftools \
view \
--output $path_file_product_bcf \
--output-type b9 \
--threads $threads \
$path_file_temporary_sort_records

# Create Tabix index for product file in VCF format with BGZip compression.
# BCFTools sometimes requires this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_product_bcf

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
