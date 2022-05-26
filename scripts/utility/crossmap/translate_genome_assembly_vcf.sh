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

# In the process of translating coordinates between genome assemblies, some SNPs
# might change chromosomes. For this reason it might be appropriate to combine
# VCF files for all chromosomes (using BCFTools "concat" command) before the
# translation and then split by chromosome after the translation (using BCFTools
# "view" command with "--regions" option).

################################################################################

################################################################################
# Organize arguments.
path_vcf_source=${1} # full path to source file in VCF format
path_vcf_product=${2} # full path to product file in VCF format
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_product_genome_assembly_sequence=${4} # full path to product genome assembly sequence file in FASTA format without compression
threads=${5} # count of processing threads to use
path_environment_crossmap=${6} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${7} # full path to installation executable of BCFTools
report=${8} # whether to print reports

###########################################################################
# Organize paths.

name_base_source_file="$(basename $path_vcf_source .vcf.gz)"
path_product_container="$(dirname $path_vcf_product)"
name_base_product_file="$(basename $path_vcf_product .vcf.gz)"
path_vcf_source_no_compression="${path_product_container}/${name_base_source_file}_no_compression.vcf"
path_vcf_product_assembly_no_compression="${path_product_container}/${name_base_product_file}_assembly_no_compression.vcf"
path_vcf_product_sort_no_compression="${path_product_container}/${name_base_product_file}_sort_no_compression.vcf"

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
--output $path_vcf_source_no_compression \
--output-type v \
--threads $threads \
$path_vcf_source

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
$path_vcf_source_no_compression \
$path_product_genome_assembly_sequence \
$path_vcf_product_assembly_no_compression

# Deactivate Virtual Environment.
deactivate
which python3

# The translation of coordinates between genome assemblies introduces records
# for SNPs in the VCF file that are not in order.
# Records for SNPs out of order in the VCF file cause a problem when trying to
# create a Tabix index.

# Read VCF file without compression in BCFTools.
# Sort coordinates in VCF file.
# Write VCF file without compression.
$path_bcftools \
sort \
--output $path_vcf_product_sort_no_compression \
--output-type v \
--temp-dir $path_product_container \
$path_vcf_product_assembly_no_compression

# Read VCF file without compression in BCFTools.
# Write VCF file with BGZip compression.
$path_bcftools \
view \
--output $path_vcf_product \
--output-type z9 \
--threads $threads \
$path_vcf_product_sort_no_compression

# Read VCF file with BGZip compression in BCFTools.
# Create Tabix index for product file in VCF format.
# BCFTools sometimes requires this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_vcf_product

# Remove temporary, intermediate files.
rm $path_vcf_source_no_compression
rm $path_vcf_product_assembly_no_compression
rm $path_vcf_product_sort_no_compression



#
