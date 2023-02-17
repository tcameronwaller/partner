#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Saevarsdottir et al, Nature, 2020 (PubMed:32581359).
# Host: https://www.decode.com/summarydata/

# Source Format
# Human Genome Assembly: GRCh37 (UK Biobank)
# Effect Allele: "A1"
# Delimiter: white space
# Columns: Chr Pos rsID A0 A1 IS-frq IS-info UKB-frq UKB-info OR-A1 P   I2  P-het
#          1   2   3    4  5  6      7       8       9        10    11  12  13    (TCW; 15 February 2023)
# Note: Meta analysis from UK and Iceland Biobanks.
# Note: Combine allele frequencies and imputation "info" scores from both cohorts.
# - - Observations in Iceland: 346,753 (cases: 4,692; controls: 342,061) of 755,406 total (45.9%)
# - - Observations in UK: 408,653 (cases: 25,542; controls: 383,111) of 755,406 total (54.1%)
# Note: combination = ((value_IS * 0.459) + (value_UKB * 0.541))

# Format Translation
# columns:

# Product Format: Team Standard
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N  Z  INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9 10 11 12   13    14

# Review: TCW; 15 Febuary 2023

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
fill_observations=${3} # logical binary indicator of whether to fill count of observations across all variants
observations=${4} # count of observations
fill_case_control=${5} # logical binary indicator of whether to fill counts of cases and controls across all variants
cases=${6} # count of cases
controls=${7} # count of controls
report=${8} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_format_${name_base_file_product}" # hopefully unique
path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_format.txt"
path_file_temporary_format_2="${path_directory_product_temporary}/${name_base_file_product}_format_2.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

# Note: TCW; 17 February 2023
# The logarithm of a negative number or zero is undefined.

# Note: TCW; 17 February 2023
# For some odd reason, the following blocks of code in awk return a syntax error
# when executed within a conditional statement.
# Comparable blocks of code in awk do not return a syntax error when executed
# without a conditional statement.
# "(a = $1); sub(/chr/, "", a); print $3, a, ..."
# "(a = $1); sub("chr", "", a); print $3, a, ..."
# "(a = $1); gsub(/chr/, "", a); print $3, a, ..."
# "(a = $1); gsub("chr", "", a); print $3, a, ..."
# Avoid the problem by removing the string in a subsequent call to awk.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" == "1" ] && [ "$fill_case_control" == "1" ]; then
  zcat $path_file_source | awk -v observations=$observations -v cases=$cases -v controls=$controls 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    if ((toupper($10) != "NA") && (($10 + 0) > 0) && (toupper($6) != "NA") && (toupper($7) != "NA") && (toupper($8) != "NA") && (toupper($9) != "NA"))
      # Calculate allele frequency and imputation quality score from non-missing values for Iceland Biobank and UK Biobank.
      print $3, $1, $2, toupper($5), toupper($4), (($6*0.459)+($8*0.541)), log($10), "NA", $11, (observations), "NA", (($7*0.459)+($9*0.541)), (cases), (controls)
    else if ((toupper($10) != "NA") && (($10 + 0) > 0) && (toupper($6) != "NA") && (toupper($7) != "NA") && (toupper($8) == "NA") && (toupper($9) == "NA"))
      # Report non-missing values from Iceland Biobank.
      print $3, $1, $2, toupper($5), toupper($4), ($6), log($10), "NA", $11, (observations), "NA", ($7), (cases), (controls)
    else if ((toupper($10) != "NA") && (($10 + 0) > 0) && (toupper($6) == "NA") && (toupper($7) == "NA") && (toupper($8) != "NA") && (toupper($9) != "NA"))
      # Report non-missing values from UK Biobank.
      print $3, $1, $2, toupper($5), toupper($4), ($8), log($10), "NA", $11, (observations), "NA", ($9), (cases), (controls)
    else if ((toupper($10) != "NA") && (($10 + 0) > 0))
      # Report missing values for allele frequency and a meaningless imputation quality score.
      print $3, $1, $2, toupper($5), toupper($4), "NA", log($10), "NA", $11, (observations), "NA", (1.0), (cases), (controls)
    else
      # Print missing values for effect, allele frequency, and a meaningless imputation quality score.
      print $3, $1, $2, toupper($5), toupper($4), "NA", "NA", "NA", $11, (observations), "NA", (1.0), (cases), (controls)
  }' >> $path_file_temporary_format
fi

# Remove "chr" prefix from chromosome identifiers.
#awk 'BEGIN {FS = " "; OFS = " "} NR > 1 { gsub(/chr/,"", $2); print } ' $path_file_temporary_format > $path_file_temporary_format
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format_2
cat $path_file_temporary_format | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  (a = $2); sub(/chr/,"", a); print $1, a, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14
} ' >> $path_file_temporary_format_2

# Compress file format.
gzip -cvf $path_file_temporary_format_2 > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after format:"
  head -10 $path_file_temporary_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
