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
path_file_vcf_source=${1} # full path to source file in VCF format
path_file_vcf_product=${2} # full path to product file in VCF format
path_file_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_file_genome_assembly_sequence=${4} # full path to genome assembly sequence file in FASTA format without compression that matches product genome assembly
threads=${5} # count of processing threads to use
path_environment_crossmap=${6} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${7} # full path to installation executable file of BCFTools
report=${8} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_vcf_product .vcf.gz)"
path_directory_product="$(dirname $path_file_vcf_product)"
path_directory_product_temporary="${path_directory_product}/temporary_mga_${name_base_file_product}" # hopefully unique

path_file_temporary_1_vcf="${path_directory_product_temporary}/${name_base_file_product}_1.vcf"
path_file_temporary_2_raw_vcf="${path_directory_product_temporary}/${name_base_file_product}_2_raw.vcf"
path_file_temporary_3_raw_bcf="${path_directory_product_temporary}/${name_base_file_product}_3_raw.bcf"

path_file_temporary_list_samples="${path_directory_product_temporary}/${name_base_file_product}_list_samples.txt"
path_file_temporary_4_sort_samples="${path_directory_product_temporary}/${name_base_file_product}_4_sort_samples.bcf"
path_file_temporary_5_sort_records="${path_directory_product_temporary}/${name_base_file_product}_5_sort_records.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

################################################################################
# Translate coordinates for chromosomes and base pair positions from human
# genome assembly GRCh38 to GRCh37.

# Read VCF file with BGZip compression in BCFTools.
# Write VCF file without compression.
$path_bcftools \
view \
--output $path_file_temporary_1_vcf \
--output-type v \
--threads $threads \
$path_file_vcf_source

# Activate Virtual Environment.
source "${path_environment_crossmap}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

# Read VCF file without compression in CrossMap.
# Translate coordinates between genome assemblies.
# Write VCF file without compression.
# CrossMap uses GZip compression ("--compress" command).
# Do not use GZip compression in order to simplify reading into BCFTools.
CrossMap.py \
vcf \
--chromid a \
$path_file_assembly_translation_chain \
$path_file_temporary_1_vcf \
$path_file_genome_assembly_sequence \
$path_file_temporary_2_raw_vcf

# Deactivate Virtual Environment.
deactivate
which python3

# Remove previous temporary file.
rm $path_file_temporary_1_vcf


# The translation of coordinates between genome assemblies introduces records
# for SNPs in the VCF file that are not in order.
# Records for SNPs out of order in the VCF file cause a problem when trying to
# create a Tabix index.

# Convert genotype files from VCF format without compression to BCF format
# without compression.
# The BCF format allows for greater performance in BCFTools.
$path_bcftools \
view \
--output $path_file_temporary_3_raw_bcf \
--output-type u \
--threads $threads \
$path_file_temporary_2_raw_vcf
# Remove previous temporary file.
rm $path_file_temporary_2_raw_vcf

# Sort samples and records in genotype files.
# BCFTools throws an error when trying to create a Tabix index if records are
# not in sort order.

# Sort samples.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf
# Extract sample identifiers from genotype file and sort these externally.
$path_bcftools \
query \
--list-samples \
$path_file_temporary_3_raw_bcf | sort > $path_file_temporary_list_samples
# Sort samples within genotype file.
$path_bcftools \
view \
--samples-file $path_file_temporary_list_samples \
--output $path_file_temporary_4_sort_samples \
--output-type u \
--threads $threads \
$path_file_temporary_3_raw_bcf
# Remove previous temporary file.
rm $path_file_temporary_3_raw_bcf


# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 4G \
--output $path_file_temporary_5_sort_records \
--output-type u \
--temp-dir $path_directory_product_temporary \
$path_file_temporary_4_sort_samples
# Remove previous temporary file.
rm $path_file_temporary_4_sort_samples

# Convert genotype files from BCF format without compression to VCF format with
# BGZip compression.
$path_bcftools \
view \
--output $path_file_vcf_product \
--output-type z9 \
--threads $threads \
$path_file_temporary_5_sort_records
# Remove previous temporary file.
rm $path_file_temporary_5_sort_records

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
rm -r $path_directory_product_temporary



#
