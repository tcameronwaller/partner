#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 23 May 2022
################################################################################
# Note

# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"

# The CrossMap tool translates coordinates for chromosomes and base pair
# positions of Single Nucleotide Polymorphisms (SNPs) from a source human genome
# assembly to a target human genome assembly.

# For translations of genotype files in Variant Call Format (VCF), CrossMap
# requires a single reference file in FASTA format for the entire target human
# genome assembly.

# Information on how to access a single FASTA file for assembly GRCh37 of the
# Human Genome: "https://www.biostars.org/p/338914/#339258"


################################################################################

# TODO: TCW; 24 May 2022
# TODO: access 'chain' file for mapping from GRCh38 to GRCh37


################################################################################
# Organize arguments.
path_dbsnp_reference=${1} # full path to file for dbSNP reference in VCF format
path_vcf_source=${2} # full path to source file in VCF format
path_vcf_product=${3} # full path to product file in VCF format
threads=${4} # count of processing threads to use
path_environment_crossmap=${5} # full path to installation of BCFTools
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
$path_translation_chain_file \
$path_vcf_source \
$path_genome_reference_fasta \
$path_vcf_product \
--chromid a \
--compress

################################################################################
# Deactivate Virtual Environment.
deactivate
which python3



#
