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

# BCFTools can read a dbSNP reference VCF file in BGZip compression format
# (BGZF) with a matching Tabix index (.tbi) for more efficient performance.
# Actually some commands in BCFTools (such as "annotate") require the matching
# Tabix index. It is possible to read a Tabix index with a non-standard name
# convention.
# Example: "bcftools view -r X:2928329 file.vcf.gz##idx##non-standard-index-name"

# It is necessary to call BCFTools separately for each command.
# For example, separate calls to BCFTools are necessary for the
# "annotate --rename-chrs" command and then the "annotate --annotations"
# command. The documentation claims that it is more efficient to pipe the
# standard output result from a first command to a second command using a symbol
# that BCFTools recognizes as standard input ("-"). In practice BCFTools seems
# to give an error reading information from the standard input because it is
# "not compressed with bgzip".
# Example:
# bcftools first-command --options source_file | bcftools second-command --options -

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
path_file_vcf_source=${1} # full path to source genotype file in VCF format
path_file_vcf_product=${2} # full path to product genotype file in VCF format
path_file_dbsnp_reference=${3} # full path to file for dbSNP reference in VCF format
threads=${4} # count of processing threads to use
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Introduce dbSNP rsID annotations VCF genotype file.
# Write to file in VCF format with BGZIP compression.

# Only introduce dbSNP rsID annotations.
$path_bcftools \
annotate \
--annotations $path_file_dbsnp_reference \
--columns ID \
--output $path_file_vcf_product \
--output-type z9 \
--threads $threads \
$path_file_vcf_source

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_vcf_product
