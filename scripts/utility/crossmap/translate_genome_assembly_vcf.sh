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

################################################################################

################################################################################
# Organize arguments.
path_vcf_source=${1} # full path to source file in VCF format
path_vcf_product=${2} # full path to product file in VCF format
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_product_genome_assembly_sequence=${4} # full path to product genome assembly sequence file in FASTA format with optional
path_environment_crossmap=${5} # full path to Python 3 environment with installation of CrossMap
report=${6} # whether to print reports

################################################################################
# Activate Virtual Environment.
source "${path_environment_crossmap}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Introduce dbSNP rsID annotations VCF genotype file.
# Write to file in VCF format with BGZIP compression.

# Translate coordinates for chromosomes and base pair positions from human
# genome assembly GRCh38 to GRCh37.

CrossMap.py \
vcf \
--chromid a \
--compress \
$path_assembly_translation_chain \
$path_vcf_source \
$path_product_genome_assembly_sequence \
$path_vcf_product

################################################################################
# Deactivate Virtual Environment.
deactivate
which python3



#
