#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# Source Format
# http://genetics.cs.ucla.edu/meta_jemdoc/
# delimiter: white space
# columns: RSID CHROM POS A1 A2 inMeta #STUDY PVALUE_FE BETA_FE STD_FE PVALUE_RE BETA_RE
# STD_RE PVALUE_RE2 STAT1_RE2 STAT2_RE2 PVALUE_BE I_SQUARE Q PVALUE_Q TAU_SQUARE
# PVALUES_OF_STUDIES(Tab_delimitered) MVALUES_OF_STUDIES(Tab_delimitered)

# Format Translation
# Use the values for fixed effects ("_FE").
# Documentation indicated that column "STD_FE" was for standard error of the
# coefficient estimate with fixed effects.
# columns: $1, $2, $3, $4, $5, "NA", $9, $10, $8, (5398), "NA", (1), "NA", "NA"

# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 23 November 2022

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
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  print $1, $2, $3, toupper($4), toupper($5), "NA", $9, $10, $8, (5398), "NA", (1), "NA", "NA"
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
