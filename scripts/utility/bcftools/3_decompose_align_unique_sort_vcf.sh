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

# Reference protocol for mapping between genomic assemblies
# "https://www.protocols.io/view/genotyping-chip-data-lift-over-to-reference-
# genome-n2bvjmbpvk5w/v2?step=1"

################################################################################

# TODO: TCW; 02 June 2022
# TODO: in addition to current quality control steps...
# TODO: 1. decompose multiallelic sites to biallelic sites (BCFTools "norm" command)
# TODO: 2. align reference alleles with the reference genome (BCFTools "norm" command, I think)

# TODO: 3. make this script versatile enough to quality control VCF files before AND after assembly mapping


################################################################################
# Organize arguments.
path_file_vcf_source=${1} # full path to source genotype file in VCF format
path_file_vcf_product=${2} # full path to product genotype file in BCF format
path_file_genome_assembly_sequence=${3} # full path to genome assembly sequence file in FASTA format without compression that matches source and product genome assembly
threads=${4} # count of processing threads to use
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_vcf_product .vcf.gz)"
path_directory_product="$(dirname $path_file_vcf_product)"
path_directory_product_temporary="${path_directory_product}/temporary_dausvcf_${name_base_file_product}" # hopefully unique

path_file_temporary_1_bcf="${path_directory_product_temporary}/${name_base_file_product}_1.bcf"
path_file_temporary_2_atom="${path_directory_product_temporary}/${name_base_file_product}_2_atom.bcf"
path_file_temporary_3_multi="${path_directory_product_temporary}/${name_base_file_product}_3_multi.bcf"
path_file_temporary_4_align="${path_directory_product_temporary}/${name_base_file_product}_4_align.bcf"
path_file_temporary_5_unique="${path_directory_product_temporary}/${name_base_file_product}_5_unique.bcf"
path_file_temporary_list_samples="${path_directory_product_temporary}/${name_base_file_product}_list_samples.txt"
path_file_temporary_6_sort_samples="${path_directory_product_temporary}/${name_base_file_product}_6_sort_samples.bcf"
path_file_temporary_7_sort_records="${path_directory_product_temporary}/${name_base_file_product}_7_sort_records.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary
rm $path_file_vcf_product

################################################################################
# Prepare genotype files for combination.

# Convert genotype files from VCF format with BGZip compression to working
# format.
# The BCF format allows for greater performance in BCFTools.
# Read file in VCF format with BGZip compression in BCFTools.
# Write file in BCF format without compression or other format.
$path_bcftools \
view \
--output $path_file_temporary_1_bcf \
--output-type u \
--threads $threads \
$path_file_vcf_source

# Decompose complex genetic features, such as multiple nucleotide variants, to
# simpler, consecutive single nucleotide variants.
$path_bcftools \
norm \
--atomize \
--atom-overlaps \* \
--output $path_file_temporary_2_atom \
--output-type u \
--threads $threads \
$path_file_temporary_1_bcf
# Remove previous temporary file.
#rm $path_file_temporary_1_bcf

# Decompose multiallelic genetic features (SNPs, etc) to separate records for
# biallelic genetic features.
$path_bcftools \
norm \
--multiallelics -any \
--output $path_file_temporary_3_multi \
--output-type u \
--threads $threads \
$path_file_temporary_2_atom
# Remove previous temporary file.
#rm $path_file_temporary_2_atom

# Align allelelic designations to the reference genome.
# bcftools norm -f ___path_to_ref_genome -c ws
# bcftools norm --fasta-ref ___path_to_ref_genome --check-ref ws
$path_bcftools \
norm \
--fasta-ref $path_file_genome_assembly_sequence \
--check-ref ws \
--output $path_file_temporary_4_align \
--output-type u \
--threads $threads \
$path_file_temporary_3_multi
# Remove previous temporary file.
#rm $path_file_temporary_3_multi

# Remove duplicate records for SNPs or other genetic features.
# I think that option "--rm-dup exact" evokes similar logic to option
# "--collapse none".
# I think both of these options require exact match of chromosome, position, and
# all alleles in order to consider records redundant.
$path_bcftools \
norm \
--rm-dup exact \
--output $path_file_temporary_5_unique \
--output-type u \
--threads $threads \
$path_file_temporary_4_align
# Remove previous temporary file.
#rm $path_file_temporary_4_align

# Sort samples.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf
# Extract sample identifiers from genotype file and sort these externally.
$path_bcftools \
query \
--list-samples \
$path_file_temporary_5_unique | sort > $path_file_temporary_list_samples
# Sort samples within genotype file.
$path_bcftools \
view \
--samples-file $path_file_temporary_list_samples \
--output $path_file_temporary_6_sort_samples \
--output-type u \
--threads $threads \
$path_file_temporary_5_unique
# Remove previous temporary file.
#rm $path_file_temporary_5_unique

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 4G \
--output $path_file_temporary_7_sort_records \
--output-type u \
--temp-dir $path_directory_product_temporary \
$path_file_temporary_6_sort_samples
# Remove previous temporary file.
#rm $path_file_temporary_6_sort_samples

# Convert genotype files from BCF format without compression to VCF format with
# BGZip compression.
$path_bcftools \
view \
--output $path_file_vcf_product \
--output-type z9 \
--threads $threads \
$path_file_temporary_7_sort_records
# Remove previous temporary file.
#rm $path_file_temporary_7_sort_records

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_vcf_product

# Remove temporary, intermediate files.
#rm -r $path_directory_product_temporary



#
