#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Pott et al, Metabolites, 2021 (PubMed:34822396).
# Host: https://zenodo.org/record/5644896

# Source Format
# Human Genome Assembly: GRCh37 (hg19)
# Effect allele: "ea"
# Delimiter: white space
# Columns: markername chr bp_hg19 ea oa eaf info nSamples nStudies beta se  p  I2  phenotype
#          1          2   3       4  5  6   7    8        9        10   11  12 13  14        (TCW; 15 February 2023)

# Format Translation
# Identifiers of SNP variants is [rsID]:[position]:[other allele]:[effect allele].
# LDSC Munge might or might not be able to interpret this original identifier format.
# Split the SNP variant identifier and extract the rs identifier.
# columns: $1, $2, $3, toupper($4), toupper($5), $6, $10, $11, $12, $8, "NA", $7, "NA", "NA"

# Product Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 23 December 2022
# check: Standard Format Columns [TCW; 22 December 2022]
# check: Study citation, PubMed, and Host website [TCW; 22 December 2022]
# check: Study field delimiters [TCW; 23 December 2022]
# check: Study source columns [TCW; 23 December 2022]
# check: Translation column order [TCW; 23 December 2022]


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

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" != "1" ] && [ "$fill_case_control" != "1" ]; then
  zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    (a = $1); split(a, b, ":"); print b[1], $2, $3, toupper($4), toupper($5), $6, $10, $11, $12, $8, "NA", $7, "NA", "NA"
  }' >> $path_file_temporary_format
fi

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

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
