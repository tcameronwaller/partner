#!/bin/bash


################################################################################
# Notes:

# National Center for Biotechnology Information (NCBI)
# Project dbSNP for Short Genetic Variations
# https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi

# The dbSNP reference from the National Center for Biotechnology Information
# (NCBI) gives information about Single Nucleotide Polymorphisms (SNPs) that
# includes unique reference SNP identifiers (rsIDs). This information is
# available in the Variant Call Format (VCF).

# Importantly, the dbSNP VCF file uses "RefSeq" identifiers for each chromosome
# in the "CHROM" column. These "RefSeq" identifiers correspond to a specific
# assembly of the human genome.
# It is convenient to translate these "RefSeq" identifiers to names of
# corresponding chromosomes.
# RefSeq
# https://www.ncbi.nlm.nih.gov/refseq/
# https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
# Select chromosomes in genome assemblies GRCh37 or GRCh38 to view their RefSeq
# identifiers in the Genome Data Browser.

# Parameter files for translation of chromosome identifiers.
# file: "translations_chromosomes_refseq_grch37p13.txt"
# file: "translations_chromosomes_refseq_grch38p14.txt"
# Review of parameter files: TCW; 6 February 2023

# Review: TCW; 6 February 2023

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
path_file_chromosomes_grch37=${2} # full path to white-space-delimited text file of translations for chromosome identifiers
path_file_chromosomes_grch38=${3} # full path to white-space-delimited text file of translations for chromosome identifiers
path_file_script=${4} # full path to file of script for translation of chromosome identifiers
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Organize paths.

path_directory_grch37_source="${path_directory_parent}/grch37"
path_directory_grch37_product="${path_directory_parent}/grch37_chromosome"
path_directory_grch38_source="${path_directory_parent}/grch38"
path_directory_grch38_product="${path_directory_parent}/grch38_chromosome"

path_file_grch37_source="${path_directory_parent}/grch37/GCF_000001405.25.gz"
path_file_grch37_product="${path_directory_parent}/grch37_chromosome/GCF_000001405.25.gz"
path_file_grch38_source="${path_directory_parent}/grch38/GCF_000001405.40.gz"
path_file_grch38_product="${path_directory_parent}/grch38_chromosome/GCF_000001405.40.gz"


# Initialize directory.
rm -r $path_directory_grch37_product
rm -r $path_directory_grch38_product
mkdir -p $path_directory_grch37_product
mkdir -p $path_directory_grch38_product

cd $path_directory_parent

###########################################################################
# Execute procedure.

# Echo each command to console.
#set -x
# Suppress echo each command to console.
#set +x

###########################################################################
# Organize directories.
# Access reference information from NCBI dbSNP.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate chromosome identifiers in dbSNP reference files for human genome assembly GRCh37."
  echo "----------"
  echo "----------"
  echo "----------"
fi
# Call script for translation.
/usr/bin/bash $path_file_script \
$path_file_grch37_source \
$path_file_grch37_product \
$path_file_chromosomes_grch37 \
4 \
$path_bcftools \
$report

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate chromosome identifiers in dbSNP reference files for human genome assembly GRCh38."
  echo "----------"
  echo "----------"
  echo "----------"
fi
# Call script for translation.
/usr/bin/bash $path_file_script \
$path_file_grch38_source \
$path_file_grch38_product \
$path_file_chromosomes_grch38 \
4 \
$path_bcftools \
$report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "translate_chromosomes_dbsnp_human_grch37_grch38_vcf.sh"
  echo "----------"
fi





#
