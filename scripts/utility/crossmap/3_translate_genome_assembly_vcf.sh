#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 23 May 2022
################################################################################
# Note

# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"

# The CrossMap tool translates coordinates for chromosomes, base-pair positions,
# and alleles of Single Nucleotide Polymorphisms (SNPs) from a source human
# genome assembly to a target human genome assembly.
# The Crossmap tool does not change annotations for SNPs (such as identifiers
# from dbSNP).

# For translations of genotype files in Variant Call Format (VCF), CrossMap
# requires a single reference file in FASTA format for the entire target human
# genome assembly.

# Information on how to access a single FASTA file for assembly GRCh37 of the
# Human Genome: "https://www.biostars.org/p/338914/#339258"

# In the process of translating coordinates between genome assemblies, some
# genetic features (SNPs etc) change chromosomes. For this reason it might be
# appropriate to combine VCF files for all chromosomes (using BCFTools "concat"
# command) before the translation and then split by chromosome after the
# translation (using BCFTools "view" command with "--regions" option).

################################################################################

################################################################################
# Organize arguments.
path_file_source_vcf=${1} # full path to source file in VCF format
path_file_product_vcf=${2} # full path to product file in VCF format
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_genome_assembly_sequence=${4} # full path to genome assembly sequence file in FASTA format without compression that matches product genome assembly
threads=${5} # count of processing threads to use
path_environment_crossmap=${6} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${7} # full path to installation executable file of BCFTools
report=${8} # whether to print reports

###########################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product_vcf .vcf.gz)"
path_directory_product="$(dirname $path_file_product_vcf)"
path_directory_product_temporary="${path_directory_product}/temporary"

path_file_temporary_decompression="${path_directory_product_temporary}/${name_base_file_product}_decompression.vcf"
path_file_temporary_assembly_raw="${path_directory_product_temporary}/${name_base_file_product}_assembly_raw.vcf"
path_file_temporary_assembly_bcf="${path_directory_product_temporary}/${name_base_file_product}_assembly_bcf.bcf"

path_file_temporary_remove_replicates="${path_directory_product_temporary}/${name_base_file_product}_replicates.bcf"
path_file_temporary_list_samples="${path_directory_product_temporary}/${name_base_file_product}_list_samples.txt"
path_file_temporary_sort_samples="${path_directory_product_temporary}/${name_base_file_product}_sort_samples.bcf"
path_file_temporary_sort_records="${path_directory_product_temporary}/${name_base_file_product}_sort_records.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Introduce dbSNP rsID annotations VCF genotype file.
# Write to file in VCF format with BGZIP compression.

# Translate coordinates for chromosomes and base pair positions from human
# genome assembly GRCh38 to GRCh37.

# CrossMap uses GZip compression.
# Do not use GZip compression in order to simplify reading into BCFTools.

# Read VCF file with BGZip compression in BCFTools.
# Write VCF file without compression.
$path_bcftools \
view \
--output $path_file_temporary_decompression \
--output-type v \
--threads $threads \
$path_file_source_vcf

# Activate Virtual Environment.
source "${path_environment_crossmap}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

# Read VCF file without compression in CrossMap.
# Translate coordinates between genome assemblies.
# Write VCF file without compression.
# Include "--compress" command to apply simple GZip compression.
CrossMap.py \
vcf \
--chromid a \
$path_assembly_translation_chain \
$path_file_temporary_decompression \
$path_genome_assembly_sequence \
$path_file_temporary_assembly_raw

# Deactivate Virtual Environment.
deactivate
which python3

# The translation of coordinates between genome assemblies introduces records
# for SNPs in the VCF file that are not in order.
# Records for SNPs out of order in the VCF file cause a problem when trying to
# create a Tabix index.

# Convert genotype files from VCF format without compression to BCF format
# without compression.
# The BCF format allows for greater performance in BCFTools.
$path_bcftools \
view \
--output $path_file_temporary_assembly_bcf \
--output-type u \
--threads $threads \
$path_file_temporary_assembly_raw

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
$path_file_temporary_assembly_bcf

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

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 4G \
--output $path_file_temporary_sort_records \
--output-type u \
--temp-dir $path_directory_product_temporary \
$path_file_temporary_sort_samples

# Convert genotype files from BCF format without compression to VCF format with
# BGZip compression.
$path_bcftools \
view \
--output $path_file_product_vcf \
--output-type z9 \
--threads $threads \
$path_file_temporary_sort_records

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_product_vcf

# Remove temporary, intermediate files.
#rm -r $path_directory_product_temporary



#
