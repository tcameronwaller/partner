#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Zhou et al, Nature Communications, 2020 (PubMed:32769997).
# Host: http://csg.sph.umich.edu/willer/public/TSH2020/
# Host: https://www.ebi.ac.uk/gwas/publications/32769997

# Source Format
# delimiter: white space
# columns: Chr Pos Allele1 Allele2 Freq N Effect SE P-value Direction HetPval

# Format Translation
# The GWAS summary statistics do not include rs identifiers for SNP variants.
# Instead, define for SNP variants a special identifier format, "chr[chromosome]_[position]_[effect allele]".
# The LDSC Munge procedure might be able to interpret this identifier.
# columns: ("chr"$1"_"$2"_"$3), $1, $2, toupper($3), toupper($4), $5, $7, $8, $9, $6, "NA", (1), "NA", "NA"

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

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ("chr"$1"_"$2"_"$3), $1, $2, toupper($3), toupper($4), $5, $7, $8, $9, $6, "NA", (1), "NA", "NA"
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
