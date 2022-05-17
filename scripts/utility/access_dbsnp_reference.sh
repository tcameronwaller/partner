#!/bin/bash


################################################################################
# Note

# National Center for Biotechnology Information (NCBI)
# Project dbSNP for Short Genetic Variations
# https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi

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

# Do not decompress the dbSNP file (BGZip format; BGZF).

################################################################################


################################################################################
# Organize arguments.
path_dbsnp_reference=${1} # full path to directory for dbSNP reference

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

###########################################################################
# Organize directories.

rm -r $path_dbsnp_reference

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_dbsnp_reference ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_dbsnp_reference
fi

cd $path_dbsnp_reference

###########################################################################
# Access reference information from NCBI dbSNP.

# Latest dbSNP build: 155
# File date: 25 May 2021
# Release: 16 June 2021
# Genome Reference Consortium (GRC) human assembly: 38 (GRCh38)

wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz           # 24 Gigabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.md5       # 512 bytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.tbi       # 3.0 Megabytes
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.tbi.md5   # 512 bytes



#
