#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Note

# BCFTools specializes in tasks on genotype files in the Variant Call Format
# (VCF) that describe Single Nucleotide Polymorphisms (SNPs) across samples.

# BCFTools is able to introduce annotation information to a target VCF file from
# a reference file also in VCF format, such as the dbSNP reference from the
# National Center for Biotechnology Information (NCBI).

# The "annotate" function of BCFTools matches records for SNPs between the two
# VCF files at least by chromosome (column: "CHROM") and position
# (column: "POS"). If the appropriate columns are present, then the function
# also matches SNPs by reference allele (column: "REF") and alternate allele
# (column: "ALT"). BCFTools requires the format of chromosome identifiers to
# match directly between both VCF files. Either both or neither VCF files must
# use the "chr" prefix in chromosome identifiers.

# It is important that both the reference and target genotype files in VCF
# format use chromosome base pair positions that correspond to the same
# Genome Reference Consortium (GRC) human assembly, such as GRCh38.

# BCFTools can read a dbSNP reference VCF file in BGZF (bgzip) compression
# format with a Tabix index (.tbi) for more efficient performance.

# BCFTools documentation
# https://samtools.github.io/bcftools/bcftools.html

# Examples of calls to the "annotate" function of BCFTools.
# http://samtools.github.io/bcftools/howtos/annotate.html

# BGZIP documentation
# http://www.htslib.org/doc/bgzip.html

# "Cheat sheet" on BCFTools
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

# Example of simple script call format for using BCFTools Annotate to introduce rsIDs
# https://gist.github.com/obenshaindw/99440ac43de07548e453
# "/usr/bin/htslib/bcftools/bcftools annotate -a /reference/dbsnp_137.b37.vcf.gz -c ID vcf_to_add_id_to.vcf"

################################################################################

################################################################################
# Organize arguments.
path_dbsnp_reference=${1} # full path to file for dbSNP reference in VCF format
path_vcf_source=${2} # full path to source file in VCF format
path_vcf_product=${3} # full path to product file in VCF format
path_bcftools=${4} # full path to installation of BCFTools
report=${5} # whether to print reports

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Introduce dbSNP rsID annotations VCF genotype file.

# VCF file format with compression
# --output-type z9 \

# Only introduce dbSNP rsID annotations.
$path_bcftools \
annotate \
--annotations $path_dbsnp_reference \
--columns ID \
--output $path_vcf_product \
--output-type v \
--threads 4 \
$path_vcf_source

# Both remove "chr" prefix from chromosome identifiers and introduce dbSNP rsID
# annotations. Only the removal of "chr" prefix works (TCW; 17 May 2022).
#$path_bcftools \
#annotate \
#--rename-chrs $path_chromosome_translations \
#--annotations $path_dbsnp_reference \
#--columns ID \
#--output $path_vcf_product \
#--output-type b9 \
#--threads 4 \
#$path_vcf_source
