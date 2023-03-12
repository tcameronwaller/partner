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

# The reference SNP cluster identifier (rsID) is a semi-stable designation of a
# genomic locus that dbSNP defines. The rsID itself is not dependent on assembly
# of the human genome (GRCh37 or GRCh38); however, for the purposes of matching
# an rsID to genomic coordinates, the dbSNP definitions do specify genomic
# coordinates according to a specific assembly (GRCh37 or GRCh38).

# Here is an informatie description of the rsID.
# Site: https://annovar.openbioinformatics.org/en/latest/articles/dbSNP/

# File Transfer Protocol (FTP)
# https://ftp.ncbi.nih.gov/snp/
# https://ftp.ncbi.nlm.nih.gov/snp/
# Information about files and formats within FTP.
# https://ftp.ncbi.nlm.nih.gov/snp/00readme.txt
# "/latest_release/"
# "Contains the most recent release of human SNP data, in VCF
# and API JSON format, along with the release notes..."
# "/latest_release/VCF/"
# "RefSNP VCF files for GRC (Genome Reference Consortium) human assembly 37
# (GCF_000001405.25) and 38 (GCF_000001405.40). Files are compressed by bgzip
# and with the tabix index."

# A copy of the dbSNP reference is also available from Human Genome Resources
# at NCBI.
# https://www.ncbi.nlm.nih.gov/genome/guide/human/#download
# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/

# Do not decompress the dbSNP file (BGZip format; BGZF).

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

# Review: TCW; 6 February 2023
# Last Accession: __ February 2023 (TCW; to NCSA server)

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_grch37="${path_directory_parent}/grch37"
path_directory_grch38="${path_directory_parent}/grch38"

# Initialize directory.
rm -r $path_directory_grch37
rm -r $path_directory_grch38
mkdir -p $path_directory_grch37
mkdir -p $path_directory_grch38

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
  echo "Accessing dbSNP reference files for human genome assembly GRCh37."
  echo "----------"
  echo "----------"
  echo "----------"
fi
cd $path_directory_grch37
# Genome Reference Consortium (GRC) human assembly: GRCh37.p13
# RefSeq accession: GCF_000001405.25
# Latest dbSNP build: 155
# Release: 16 June 2021
# File date: 16 November 2022
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz           # 24 Gigabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.md5       # 54 bytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi       # 2.9 Megabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi.md5   # 58 bytes

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Accessing dbSNP reference files for human genome assembly GRCh38."
  echo "----------"
  echo "----------"
  echo "----------"
fi
cd $path_directory_grch38
# Genome Reference Consortium (GRC) human assembly: GRCh38.p14
# RefSeq accession: GCF_000001405.40
# Latest dbSNP build: 155
# Release: 16 June 2021
# File date: 16 November 2022
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz           # 25 Gigabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.md5       # 54 bytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi       # 3.0 Megabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi.md5   # 58 bytes

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "access_dbsnp_human_grch37_grch38_vcf.sh"
  echo "----------"
fi





#
