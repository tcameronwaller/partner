#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Teumer et al, Nature Communications, 2018 (PubMed:30367059).
# Host: https://transfer.sysepi.medizin.uni-greifswald.de/thyroidomics/datasets/

# Source Format
# Human Genome Assembly: GRCh37 (TCW; 24 January 2023)
# effect allele: "Allele1" ("coding allele") (TCW; 24 January 2023)
# delimiter: comma
# columns: " MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,P.value,N,I2 " # (TCW; 24 January 2023)
# Note:
# The source format is the same for linear and logistic GWAS in collection (TCW; 24 January 2023).
# Format of "MarkerName" is "[chromosome]:[position]:[type such as SNP]" (TCW; 24 January 2023).

# Format Translation
# Identifiers of SNP variants is [chromosome]:[position]:"SNP".
# Extract the chromosome and position.
# Define a temporary identifier [chromosome]:[position]:"SNP".

# TODO: TCW; 20 January 2023

# columns: $1, [extract], [extract], toupper($2), toupper($3), $4, $5, $6, $7, $8, "NA", (1), "NA", "NA"

# Product Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 24 January 2023


###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

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

# simple: print $1, $2, $3, toupper($4), toupper($5), $6, $10, $11, $12, $8, "NA", $7, "NA", "NA"

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
zcat $path_file_source | awk 'BEGIN {FS = ","; OFS = " "} NR > 1 {
  (a = $1); split(a, b, ":"); print a, b[1], b[2], toupper($2), toupper($3), $4, $5, $6, $7, $8, "NA", (1.0), "NA", "NA"
}' >> $path_file_temporary_format

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
