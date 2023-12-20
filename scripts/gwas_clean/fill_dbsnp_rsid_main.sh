#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 December 2023
# Date, last execution: 19 December 2023
# Review: 19 December 2023
################################################################################
# Notes:

# This script extracts reference SNP cluster identifiers (rsIDs) from the dbSNP
# reference and matches them to genomic coordinates in GWAS summary statistics
# to fill missing rsIDs in the GWAS summary statistics.

# Allele designations.
# 1. In the dbSNP VCF file, the reference allele ("REF") corresponds to the
# reference human genome, and the alternate allele(s) ("ALT") correspond to the
# alternative allele(s) that sometimes replaces the reference allele at that
# position in the human genome.
# 2. In GWAS summary statistics, either the dbSNP reference allele or the dbSNP
# alternate allele might correspond to the effect allele. Whichever allele is
# the effect allele corresponds to the directionality of the GWAS association
# effect parameter (beta).

# In situations for which dbSNP references multiple rsIDs for a single SNP, this
# procedure keeps only the first rsID from the list.

# The script below is a companion to perform the first few operations of this
# script to save time in subsequent iterations.
# "/.../partner/scripts/bcftools/extract_dbsnp_biallelic_sites_allele_identifiers.sh"

# TODO: TCW; 13 December 2023
# Write out two different "product files"...
# 1. Filtered to include only SNPs that had matching rsIDs from dbSNP
#   - This table is useful for a count of the "successful" SNPs if nothing else.
# 2. All SNPs from the original GWAS sum stats regardless.
#   - Subsequent procedures might be able to use more SNPs.
# IDEA: --> introduce an argument "strict" to determine whether to write out all SNPs or only matching


##########
# Standard format of GWAS summary statistics.
# This is the obligatory format of the source and product GWAS summary
# statistics.
# File name suffix: ".txt.gz"
# File compression: Gzip
# Delimiter: white space (" ")
# Columns and sequence:
# $1  $2  $3 $4 $5 $6   $7   $8 $9 $10 $11 $12  $13   $14
# SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT

##########
# Limitations.
# 1. The current implementation of this script makes internal reference to
# anonymous hard-coded file paths for the sake of convenience when executing in
# a specific computational environment. The script could of course be adapted to
# accept these file path variables as parameters (arguments).
# 2. For simplicity, this procedure filters to biallelic sites before extracting
# information about SNPs and their rsIDs from the dbSNP reference.

##########
# Testing.
# TCW; 13 December 2023
# Version of dbSNP: build 155
# Accession of dbSNP: TCW; 2023-02-06
# Translation of chromosome designations in dbSNP: TCW; 2023-02-15
# I tested the procedure using the GWAS summary statistics for levels of
# testosterone in females from Ruth et al, 2020 (PubMed:32042192). This set of
# GWAS summary statistics already included rsIDs for most SNPs, and 99.765% of
# the rsIDs matched from dbSNP were identical to the original rsIDs. Hence this
# procedure has high accuracy in its matching of rsIDs from dbSNP. The main
# limitation of this procedure is the accommodation of biallelic sites only from
# dbSNP, and the test GWAS summary statistics had lost about half (47%) of their
# SNPs after this procedure. So in summary, this procedure matches rsIDs from
# dbSNP with high accuracy but also loses a considerable proportion of SNPs,
# many of which might be for multiallelic sites.


################################################################################


################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
strict=${3} # whether to return GWAS summary statistics filtered to SNPs with successful match to dbSNP rsID
report=${4} # whether to print reports


################################################################################
# Organize paths.

# Tools.
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_dbsnp="${path_directory_reference}/dbsnp/grch37_chromosome" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_product="$(dirname $path_file_gwas_product)"
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
path_directory_temporary="${path_directory_product}/temporary_2713956_${name_base_file_gwas_product}" # must be unique

# Files.
path_file_reference_dbsnp="${path_directory_reference_dbsnp}/GCF_000001405.25.gz" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15

# Scripts.
path_file_script_1="${path_directory_partner_scripts}/gwas_clean/fill_dbsnp_rsid_1.sh"
path_file_script_2="${path_directory_partner_scripts}/gwas_clean/fill_dbsnp_rsid_2.sh"
path_file_script_3="${path_directory_partner_scripts}/gwas_clean/fill_dbsnp_rsid_3.sh"
path_file_script_4="${path_directory_partner_scripts}/gwas_clean/fill_dbsnp_rsid_4.sh"

# Initialize files.
rm $path_file_gwas_product

# Initialize directories.
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product

###########################################################################
# Execute procedure.

##########
# Call script 1.
# Filter dbSNP to biallelic sites.
# Export dbSNP to flat text table.
# Organize the dbSNP tables for merge 1 and merge 2.
if false; then
  /usr/bin/bash $path_file_script_1 \
  $path_file_reference_dbsnp \
  $path_directory_temporary \
  $path_bcftools \
  $report
fi

##########
# Call script 2.
# Merge 1 and subsequent clean up.
/usr/bin/bash $path_file_script_2 \
$path_file_gwas_source \
$path_directory_reference_dbsnp \
$path_directory_temporary \
$report

##########
# Call script 3.
# Merge 2 and subsequent clean up.
/usr/bin/bash $path_file_script_3 \
$path_directory_temporary \
$report

##########
# Call script 4.
# Integration, format translation, and writing product file.
/usr/bin/bash $path_file_script_4 \
$path_file_gwas_product \
$path_directory_temporary \
$strict \
$report

##########
# Call script 5.
# Reports.
/usr/bin/bash $path_file_script_5 \
$path_file_gwas_source \
$path_file_gwas_product \
$path_directory_temporary \
$strict \
$report



################################################################################
# Remove temporary, intermediate directories and files.

# Suppress this block for debugging.
if true; then
  rm -r $path_directory_temporary
fi



#
