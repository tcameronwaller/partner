#!/bin/bash

################################################################################
# Note


# The "annotate" function of BCFTools matches records for SNPs between the two
# VCF files at least by chromosome (column: "CHROM") and position
# (column: "POS").
# If the appropriate columns are present, then the function also matches by
# reference (column: "REF") and alternate (column: "ALT") alleles.


# Examples of calls to the "annotate" function of BCFTools.
# http://samtools.github.io/bcftools/howtos/annotate.html


# Specify desired version of dbSNP reference using cooridnates from appropriate
# assembly of the human genome.

# BCFTools can work with dbSNP reference vcf file in BGZF (bgzip) format with a
# Tabix index (.tbi) for more efficient performance.


# BGZIP documentation
# http://www.htslib.org/doc/bgzip.html

# BCFTools documentation
# https://samtools.github.io/bcftools/bcftools.html

# Example of simple script call format for using BCFTools Annotate to introduce rsIDs
# https://gist.github.com/obenshaindw/99440ac43de07548e453
# "/usr/bin/htslib/bcftools/bcftools annotate -a /reference/dbsnp_137.b37.vcf.gz -c ID vcf_to_add_id_to.vcf"


# "Cheat sheet" on BCFTools
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

# "bcftools annotate -c 'ID' -a source.vcf.gz target.vcf.gz"

################################################################################

################################################################################
# Organize arguments.
path_dbsnp_reference=${1} # full path to directory for dbSNP reference
path_vcf_source=${2} # full path to source file in VCF format
path_vcf_product=${3} # full path to product file in VCF format
path_bcftools=${4} # full path to installation of BCFTools
report=${5} # whether to print reports

################################################################################
# Introduce annotation information from dbSNP reference to file in VCF format.

$path_bcftools \
annotate \
--annotations $path_dbsnp_reference \
--columns ID $path_vcf_source \
--output $path_vcf_product \
--output-type b9 \
--threads 4
