#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Pott et al, Journal of Clinical Endocrinology and Metabolism, 2019 (PubMed:31169883).
# Host: https://www.health-atlas.de/studies/19
# Host: https://www.health-atlas.de/data_files/207

# Source Format
# delimiter: white space
# columns ("GWAS"): markername chr bp_hg19 effect_allele other_allele effect_allele_freq info n beta se p
# columns ("META"): markername chr bp_hg19 effect_allele other_allele effect_allele_freq min_info n beta se p CochransQ pCochransQ

# Format Translation
# Identifiers of SNP variants is [rsID]:[position]:[other allele]:[effect allele].
# LDSC Munge might or might not be able to interpret this original identifier format.
# Split the SNP variant identifier and extract the rs identifier.
# columns: $1, $2, $3, toupper($4), toupper($5), $6, $9, $10, $11, $8, "NA", $7, "NA", "NA"

# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; __ December 2022

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

# simple: print $1, $2, $3, toupper($4), toupper($5), $6, $9, $10, $11, $8, "NA", $7, "NA", "NA"

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  (a = $1); split(a, b, ":"); print b[1], $2, $3, toupper($4), toupper($5), $6, $9, $10, $11, $8, "NA", $7, "NA", "NA"
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
