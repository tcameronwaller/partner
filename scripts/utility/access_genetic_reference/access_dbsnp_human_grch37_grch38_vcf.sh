#!/bin/bash


################################################################################
# Note

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
# (GCF_000001405.25) and 38 (GCF_000001405.39). Files are compressed by bgzip
# and with the tabix index."

# A copy of the dbSNP reference is also available from Human Genome Resources
# at NCBI.
# https://www.ncbi.nlm.nih.gov/genome/guide/human/#download
# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/

# Do not decompress the dbSNP file (BGZip format; BGZF).

################################################################################

################################################################################
# Organize arguments.
human_genome_assembly=${1} # human genome assembly, either 'grch37' or 'grch38'
path_dbsnp_parent_container=${2} # full path to parent directory for dbSNP references

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

###########################################################################
# Organize directories.
# Access reference information from NCBI dbSNP.

if [[ "$human_genome_assembly" == "grch37" ]]; then
  echo "Accessing dbSNP reference files for human genome assembly GRCh37."
  path_dbsnp_container="${path_dbsnp_parent_container}/grch37_raw"
  rm -r $path_dbsnp_container
  mkdir -p "${path_dbsnp_container}"
  cd $path_dbsnp_container
  # Genome Reference Consortium (GRC) human assembly: GRCh37.p13
  # RefSeq accession: GCF_000001405.25
  # Latest dbSNP build: 155
  # File date: 25 May 2021
  # Release: 16 June 2021
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz           # 23 Gigabytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.md5       # 512 bytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi       # 2.9 Megabytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi.md5   # 512 bytes
elif [[ "$human_genome_assembly" == "grch38" ]]; then
  echo "Accessing dbSNP reference files for human genome assembly GRCh38."
  path_dbsnp_container="${path_dbsnp_parent_container}/grch38_raw"
  rm -r $path_dbsnp_container
  mkdir -p "${path_dbsnp_container}"
  cd $path_dbsnp_container
  # Genome Reference Consortium (GRC) human assembly: GRCh38.p14
  # RefSeq accession: GCF_000001405.39
  # Latest dbSNP build: 155
  # File date: 25 May 2021
  # Release: 16 June 2021
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz           # 24 Gigabytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.md5       # 512 bytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.tbi       # 3.0 Megabytes
  wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.tbi.md5   # 512 bytes
else
  echo "invalid specification of human genome assembly"
fi



#
